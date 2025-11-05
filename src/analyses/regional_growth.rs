use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use ml_helpers::other_regression::power_law_with_intercept::{solve, PowerLawIntercept};

use crate::{
    analyses::regional_helpers::get_edge_windows,
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, Orientation, PathSegment},
    html_report::{AnalysisSection, ReportItem, Window},
    util::{get_default_plot_downloads, CountType},
};

use super::{
    regional_helpers::{get_windows, split_ref_paths},
    Analysis, ConstructibleAnalysis, InputRequirement,
};

pub struct RegionalGrowth {
    reference_text: PathSegment,
    window_size: usize,
    count_type: CountType,
    values: HashMap<PathSegment, Vec<Window>>,
    coverage: usize,
    log_windows: bool,
}

impl Analysis for RegionalGrowth {
    fn generate_table(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        let gb = gb.expect("Regional growth should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let mut text = format!("track type=bedGraph name=\"Panacus growth\" description=\"Average growth for a window of bps\" visibility=full color=200,100,0 priority=20\n");
        for (sequence_id, sequence) in &self.values {
            for window in sequence {
                let line = format!(
                    "{} {} {} {}\n",
                    sequence_id.to_string(),
                    window.start,
                    window.end,
                    window.values[0]
                );
                text.push_str(&line);
            }
        }
        Ok(text)
    }

    fn generate_report_section(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let gb = gb.expect("Regional growth should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let id_prefix = format!(
            "regional-growth-{}",
            self.get_run_id(gb)
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let all_values = self.values.clone();
        let items: Vec<ReportItem> = all_values
            .into_iter()
            .map(|(sequence, values)| ReportItem::Chromosomal {
                id: format!("{id_prefix}-{}-{}", self.count_type, sequence.to_string()),
                name: gb.get_fname(),
                labels: vec![
                    "Growth".to_string(),
                    "c".to_string(),
                    "Adapted c".to_string(),
                    "Successfull fit".to_string(),
                ],
                is_diverging: true,
                contains_outliers: false,
                sequence: sequence.to_string(),
                values,
            })
            .collect();
        let table_text = self.generate_table(Some(gb))?;
        let table_text = format!("`{}`", table_text);
        let regional_growth_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{}", self.count_type),
            analysis: "Regional Growth".to_string(),
            table: Some(table_text),
            run_name: self.get_run_name(gb),
            run_id: self.get_run_id(gb),
            countable: self.count_type.to_string(),
            items,
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(regional_growth_tabs)
    }

    fn get_graph_requirements(&self) -> std::collections::HashSet<super::InputRequirement> {
        let mut req = HashSet::from([
            InputRequirement::Hist,
            InputRequirement::Bp,
            InputRequirement::Node,
            InputRequirement::Path(self.reference_text.clone()),
        ]);
        req.extend(InputRequirement::from_count(self.count_type).into_iter());
        req
    }

    fn get_type(&self) -> String {
        "Regionalgrowth".to_string()
    }
}

impl ConstructibleAnalysis for RegionalGrowth {
    fn from_parameter(parameter: crate::analysis_parameter::AnalysisParameter) -> Self {
        let (reference, window_size, count_type, coverage, log_windows) = match parameter {
            AnalysisParameter::RegionalGrowth {
                reference,
                window_size,
                count_type,
                coverage,
                log_windows,
            } => (reference, window_size, count_type, coverage, log_windows),
            _ => panic!("Regional growth should only be called with correct parameter"),
        };
        Self {
            reference_text: PathSegment::from_str(&reference),
            count_type,
            window_size,
            values: HashMap::new(),
            coverage,
            log_windows,
        }
    }
}

impl RegionalGrowth {
    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-regionalgrowth", gb.get_run_id())
    }

    fn set_values(&mut self, gb: &GraphBroker) {
        match self.count_type {
            CountType::Node | CountType::Bp => self.set_values_nodes(gb),
            CountType::Edge => self.set_values_edges(gb),
            CountType::All => {
                unimplemented!("Regional growth for count type all has not been implemented yet")
            }
        }
        let sum_lens = self
            .values
            .iter()
            .map(|(_s, v)| v.iter().map(|w| (w.end - w.start) as u64).sum::<u64>())
            .sum::<u64>() as f64;
        let sum_c = self
            .values
            .iter()
            .map(|(_s, v)| {
                v.iter()
                    .map(|w| (w.end - w.start) as f64 * w.values[1])
                    .sum::<f64>()
            })
            .sum::<f64>();
        let sum_alpha = self
            .values
            .iter()
            .map(|(_s, v)| {
                v.iter()
                    .map(|w| (w.end - w.start) as f64 * w.values[0])
                    .sum::<f64>()
            })
            .sum::<f64>();
        let avg_c = sum_c / sum_lens;
        let avg_alpha = sum_alpha / sum_lens;
        log::info!(
            "Calculated growth for {} with cov = {}, avg r2 = {}, avg alpha = {}",
            self.count_type,
            self.coverage,
            avg_c,
            avg_alpha
        );
    }

    fn set_values_nodes(&mut self, gb: &GraphBroker) {
        let edge2id = gb.get_edges();
        let neighbors: HashMap<(ItemId, Orientation), HashSet<ItemId>> = edge2id
            .keys()
            .flat_map(|x| [((x.0, x.1), x.2), ((x.2, x.3), x.0)])
            .fold(HashMap::new(), |mut acc, (k, v)| {
                acc.entry(k).or_default().insert(v);
                acc
            });
        let ref_paths = gb.get_all_matchings_paths(&self.reference_text);
        let ref_paths = split_ref_paths(ref_paths);
        let node_lens = gb.get_node_lens();
        let mut all_growths_of_windows: HashMap<PathSegment, Vec<Window>> = HashMap::new();
        for (sequence_id, sequence) in ref_paths {
            for (contig_id, contig) in sequence {
                let contig_start = contig_id.start.unwrap_or_default();
                let windows = get_windows(
                    contig,
                    node_lens,
                    self.window_size,
                    &neighbors,
                    contig_start,
                    self.log_windows,
                );
                log::info!("Calculating growth for {} windows", windows.len());
                let growths_of_windows: Vec<Window> = windows
                    // .par_iter()
                    .iter()
                    .filter_map(|(window, start, end)| {
                        let params: Vec<f64> = {
                            let indices: HashSet<usize> =
                                window.iter().map(|(idx, _)| idx.0 as usize).collect(); // Necessary to first create HashSet to remove duplicate nodes
                            let indices: Vec<usize> = indices.into_iter().collect();
                            if indices.len() == 1 {
                                return None;
                            }
                            let uncovered_bps = if self.count_type == CountType::Bp {
                                let uncovered_bps = window
                                    .iter()
                                    .filter_map(|(node, length)| {
                                        if *length == 0
                                            || *length == node_lens[node.0 as usize] as usize
                                        {
                                            None
                                        } else {
                                            Some((
                                                node.0,
                                                node_lens[node.0 as usize] as usize - *length,
                                            ))
                                        }
                                    })
                                    .collect::<HashMap<u64, usize>>();
                                Some(uncovered_bps)
                            } else {
                                None
                            };
                            let hist = gb.get_hist_for_subset(
                                self.count_type,
                                &indices,
                                uncovered_bps.clone(),
                            );
                            let x: Vec<f64> = (1..hist.len()).map(|x| (x as f64)).collect();
                            let y = hist.into_iter().skip(1).collect::<Vec<f64>>();
                            let growth = gb.get_growth_for_subset(
                                self.count_type,
                                &indices,
                                uncovered_bps,
                                self.coverage,
                            );
                            let f_new: Vec<f64> =
                                growth.iter().tuple_windows().map(|(a, b)| b - a).collect();
                            let f_new_x: Vec<f64> =
                                (2..f_new.len() + 2).map(|x| x as f64).collect();
                            if self.log_windows {
                                let node_ids_to_names = gb.get_node_names();
                                let nodes_text = window
                                    .iter()
                                    .map(|(idx, _)| *idx)
                                    .collect::<HashSet<ItemId>>()
                                    .into_iter()
                                    .map(|id| str::from_utf8(&node_ids_to_names[&id]).unwrap())
                                    .join("\t");
                                eprintln!(
                                    "N {}-{}\t{}",
                                    start + contig_start,
                                    end + contig_start,
                                    nodes_text
                                );
                                x.iter().zip(y.iter()).zip(growth.iter()).for_each(
                                    |((x_val, y_val), growth_val)| {
                                        eprintln!(
                                            "W {}-{}\t{}\t{}\t{}",
                                            start + contig_start,
                                            end + contig_start,
                                            x_val,
                                            y_val,
                                            growth_val,
                                        );
                                    },
                                );
                            }
                            let power_law =
                                PowerLawIntercept::from(f_new_x[1..].to_vec(), f_new[1..].to_vec());
                            let params = solve(power_law);
                            let alpha = params[1];
                            let c = params[2];
                            let c_adapted = c / params[0];
                            let success_bit = params[3];
                            if success_bit < 0.5 {
                                return None;
                            }
                            vec![alpha, c, c_adapted, success_bit]
                        };
                        let window = Window {
                            start: *start + contig_start,
                            end: *end + contig_start,
                            values: params,
                        };
                        Some(window)
                    })
                    .collect();
                all_growths_of_windows
                    .entry(sequence_id.clone())
                    .or_default()
                    .extend(growths_of_windows);
            }
        }
        self.values = all_growths_of_windows;
    }

    fn set_values_edges(&mut self, gb: &GraphBroker) {
        let edge2id = gb.get_edges();
        let neighbors: HashMap<(ItemId, Orientation), HashSet<(ItemId, Orientation)>> = edge2id
            .keys()
            .flat_map(|x| [((x.0, x.1), (x.2, x.3)), ((x.2, x.3), (x.0, x.1))])
            .fold(HashMap::new(), |mut acc, (k, v)| {
                acc.entry(k).or_default().insert(v);
                acc
            });
        let ref_paths = gb.get_all_matchings_paths(&self.reference_text);
        let ref_paths = split_ref_paths(ref_paths);
        let node_lens = gb.get_node_lens();
        let mut all_growths_of_windows: HashMap<PathSegment, Vec<Window>> = HashMap::new();
        for (sequence_id, sequence) in ref_paths {
            for (contig_id, contig) in sequence {
                let contig_start = contig_id.start.unwrap_or_default();
                let windows = get_edge_windows(
                    contig,
                    node_lens,
                    self.window_size,
                    &neighbors,
                    &edge2id,
                    contig_start,
                );
                log::info!("Calculating growth for {} windows", windows.len());
                let growths_of_windows: Vec<Window> = windows
                    .par_iter()
                    .map(|(window, start, end)| {
                        let params = {
                            let indeces: Vec<usize> =
                                window.iter().map(|idx| idx.0 as usize).collect();
                            let uncovered_bps = None;
                            let growth = gb.get_growth_for_subset(
                                self.count_type,
                                &indeces,
                                uncovered_bps,
                                self.coverage,
                            );
                            let x: Vec<f64> = (1..=growth.len())
                                .map(|x| (x as f64).log10())
                                .map(|x| {
                                    if x.is_infinite() && x.is_sign_negative() {
                                        -1000.0
                                    } else {
                                        x
                                    }
                                })
                                .collect();
                            let y: Vec<f64> = vec![0.0]
                                .into_iter()
                                .chain(growth.into_iter())
                                .tuple_windows()
                                .map(|(x_prev, x_curr)| (x_curr - x_prev).log10())
                                .map(|x| {
                                    if x.is_infinite() && x.is_sign_negative() {
                                        -1000.0
                                    } else {
                                        x
                                    }
                                })
                                .collect();
                            let huber = PowerLawIntercept::from(x.clone(), y.clone());
                            let params = solve(huber);
                            let r_squared = calculate_r_squared(&x, &y, params[0], params[1]);
                            vec![-params[0], r_squared]
                        };
                        let window = Window {
                            start: *start + contig_start,
                            end: *end + contig_start,
                            values: params,
                        };
                        window
                    })
                    .collect();
                all_growths_of_windows
                    .entry(sequence_id.clone())
                    .or_default()
                    .extend(growths_of_windows);
            }
        }
        self.values = all_growths_of_windows;
    }
}

pub fn calculate_r_squared(x: &Vec<f64>, y: &Vec<f64>, m: f64, b: f64) -> f64 {
    let f: Vec<f64> = x.iter().map(|x_val| m * x_val + b).collect();
    let ss_res: f64 = y
        .iter()
        .zip(f.iter())
        .map(|(y_val, f_val)| (y_val - f_val).powi(2))
        .sum();
    let y_hat: f64 = y.iter().copied().sum::<f64>() / y.len() as f64;
    let ss_tot: f64 = y.iter().map(|y_val| (y_val - y_hat).powi(2)).sum();
    1.0 - (ss_res / ss_tot)
}
