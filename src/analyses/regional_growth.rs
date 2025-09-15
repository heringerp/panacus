use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use ml_helpers::linear_regression::huber_regressor::{solve, HuberRegressor};

use crate::{
    analyses::regional_helpers::get_edge_windows,
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, Orientation, PathSegment},
    html_report::{AnalysisSection, ReportItem},
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
    values: HashMap<PathSegment, Vec<(f64, usize, usize)>>,
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
            for (growth, start, end) in sequence {
                let line = format!("{} {} {} {}\n", sequence_id.to_string(), start, end, growth,);
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
                label: "Average growth".to_string(),
                is_diverging: true,
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
        let node_names = if self.log_windows {
            gb.get_node_names()
        } else {
            HashMap::new()
        };
        let mut all_growths_of_windows: HashMap<PathSegment, Vec<(f64, usize, usize)>> =
            HashMap::new();
        for (sequence_id, sequence) in ref_paths {
            for (contig_id, contig) in sequence {
                let contig_start = contig_id.start.unwrap_or_default();
                let windows = get_windows(
                    contig,
                    node_lens,
                    self.window_size,
                    &neighbors,
                    contig_start,
                );
                log::info!("Calculating growth for {} windows", windows.len());
                let growths_of_windows: Vec<(f64, usize, usize)> = windows
                    // .par_iter()
                    .iter()
                    .filter_map(|(window, start, end)| {
                        Some((
                            {
                                if self.log_windows {
                                    eprint!(
                                        "window {}-{}\t",
                                        start + contig_start,
                                        end + contig_start
                                    );
                                    let text = window
                                        .iter()
                                        .map(|(id, _)| {
                                            str::from_utf8(&node_names[id])
                                                .expect("Node name can be parsed to string")
                                        })
                                        .join("\t");
                                    eprintln!("{}", text);
                                }
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
                                let growth = gb.get_growth_for_subset(
                                    self.count_type,
                                    &indices,
                                    uncovered_bps,
                                    self.coverage,
                                );
                                let x: Vec<f64> = (1..=growth.len())
                                    .map(|x| (x as f64).log10())
                                    .map(|x| {
                                        if x.is_infinite() && x.is_sign_negative() {
                                            -100000.0
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
                                            -100000.0
                                        } else {
                                            x
                                        }
                                    })
                                    .collect();
                                let huber = HuberRegressor::from(x, y);
                                let params = solve(huber);
                                -params[0]
                            },
                            *start + contig_start,
                            *end + contig_start,
                        ))
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
        let mut all_growths_of_windows: HashMap<PathSegment, Vec<(f64, usize, usize)>> =
            HashMap::new();
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
                let growths_of_windows: Vec<(f64, usize, usize)> = windows
                    .par_iter()
                    .map(|(window, start, end)| {
                        (
                            {
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
                                let huber = HuberRegressor::from(x, y);
                                let params = solve(huber);
                                -params[0]
                            },
                            *start + contig_start,
                            *end + contig_start,
                        )
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
