use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

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

pub struct RegionalCount {
    reference_text: PathSegment,
    window_size: usize,
    count_type: CountType,
    values: HashMap<PathSegment, Vec<(f64, usize, usize)>>,
}

impl Analysis for RegionalCount {
    fn generate_table(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        let gb = gb.expect("Regional count should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let mut text = format!("track type=bedGraph name=\"Panacus count\" description=\"Average count for a window of bps\" visibility=full color=200,100,0 priority=20\n");
        for (sequence_id, sequence) in &self.values {
            for (count, start, end) in sequence {
                let line = format!("{} {} {} {}\n", sequence_id.to_string(), start, end, count);
                text.push_str(&line);
            }
        }
        Ok(text)
    }

    fn generate_report_section(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let gb = gb.expect("Regional count should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let id_prefix = format!(
            "regional-count-{}",
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
                label: "Count".to_string(),
                second_label: "".to_string(),
                is_diverging: false,
                sequence: sequence.to_string(),
                values: values.into_iter().map(|(v, s, e)| (v, 0.0, s, e)).collect(),
            })
            .collect();
        let table_text = self.generate_table(Some(gb))?;
        let table_text = format!("`{}`", table_text);
        let regional_count_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{}", self.count_type),
            analysis: "Regional Count".to_string(),
            table: Some(table_text),
            run_name: self.get_run_name(gb),
            run_id: self.get_run_id(gb),
            countable: self.count_type.to_string(),
            items,
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(regional_count_tabs)
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
        "RegionalCount".to_string()
    }
}

impl ConstructibleAnalysis for RegionalCount {
    fn from_parameter(parameter: crate::analysis_parameter::AnalysisParameter) -> Self {
        let (reference, window_size, count_type) = match parameter {
            AnalysisParameter::RegionalCount {
                reference,
                window_size,
                count_type,
            } => (reference, window_size, count_type),
            _ => panic!("Regional count should only be called with correct parameter"),
        };
        Self {
            reference_text: PathSegment::from_str(&reference),
            count_type,
            window_size,
            values: HashMap::new(),
        }
    }
}

impl RegionalCount {
    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-regionalcount", gb.get_run_id())
    }

    fn set_values(&mut self, gb: &GraphBroker) {
        match self.count_type {
            CountType::Node | CountType::Bp => self.set_values_nodes(gb),
            CountType::Edge => self.set_values_edges(gb),
            CountType::All => {
                unimplemented!("Regional count for count type all has not been implemented yet")
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
        let mut all_counts_of_windows: HashMap<PathSegment, Vec<(f64, usize, usize)>> =
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
                log::info!("Calculating count for {} windows", windows.len());
                let counts_of_windows: Vec<(f64, usize, usize)> = windows
                    .par_iter()
                    .map(|(window, start, end)| {
                        (
                            {
                                let indeces: Vec<usize> =
                                    window.iter().map(|(idx, _)| idx.0 as usize).collect();
                                let bp_length = if self.count_type == CountType::Bp {
                                    let bp_length = window
                                        .iter()
                                        .map(|(node, length)| {
                                            if *length == 0 {
                                                node_lens[node.0 as usize] as usize
                                            } else {
                                                *length
                                            }
                                        })
                                        .sum();
                                    Some(bp_length)
                                } else {
                                    None
                                };
                                let count = get_count(self.count_type, &indeces, bp_length);
                                count as f64
                            },
                            *start + contig_start,
                            *end + contig_start,
                        )
                    })
                    .collect();
                all_counts_of_windows
                    .entry(sequence_id.clone())
                    .or_default()
                    .extend(counts_of_windows);
            }
        }
        self.values = all_counts_of_windows;
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
        let mut all_counts_of_windows: HashMap<PathSegment, Vec<(f64, usize, usize)>> =
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
                log::info!("Calculating count for {} windows", windows.len());
                let counts_of_windows: Vec<(f64, usize, usize)> = windows
                    .par_iter()
                    .map(|(window, start, end)| {
                        (
                            {
                                let indeces: Vec<usize> =
                                    window.iter().map(|idx| idx.0 as usize).collect();
                                let bp_length = None;
                                let count = get_count(self.count_type, &indeces, bp_length);
                                count as f64
                            },
                            *start + contig_start,
                            *end + contig_start,
                        )
                    })
                    .collect();
                all_counts_of_windows
                    .entry(sequence_id.clone())
                    .or_default()
                    .extend(counts_of_windows);
            }
        }
        self.values = all_counts_of_windows;
    }
}

fn get_count(count: CountType, indices: &Vec<usize>, bp_length: Option<usize>) -> usize {
    match count {
        CountType::Node | CountType::Edge => indices.len(),
        CountType::Bp => bp_length.expect(""),
        CountType::All => unimplemented!("Count type All not yet implemented for regional count"),
    }
}
