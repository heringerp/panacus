use std::collections::{HashMap, HashSet};

use crate::{
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, Orientation, PathSegment},
    html_report::{AnalysisSection, ReportItem},
    util::{get_default_plot_downloads, CountType},
};

use super::{
    regional_helpers::{get_windows, split_ref_paths},
    Analysis, ConstructibleAnalysis, InputRequirement,
};

pub struct RegionalDegree {
    reference_text: PathSegment,
    window_size: usize,
    values: HashMap<PathSegment, Vec<(f64, usize, usize)>>,
}

impl Analysis for RegionalDegree {
    fn generate_table(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        let gb = gb.expect("Regional degree should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let mut text = format!("track type=bedGraph name=\"Panacus Degree\" description=\"Average degree for a window of bps\" visibility=full color=200,100,0 priority=20\n");
        for (sequence_id, sequence) in &self.values {
            for (degree, start, end) in sequence {
                let line = format!("{} {} {} {}\n", sequence_id.to_string(), start, end, degree,);
                text.push_str(&line);
            }
        }
        Ok(text)
    }

    fn generate_report_section(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let gb = gb.expect("Regional degree should always have a graph");
        if self.values.is_empty() {
            self.set_values(gb);
        }
        let id_prefix = format!(
            "regional-degree-{}",
            self.get_run_id(gb)
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let all_values = self.values.clone();
        let items: Vec<ReportItem> = all_values
            .into_iter()
            .map(|(sequence, values)| ReportItem::Chromosomal {
                id: format!("{id_prefix}-{}-{}", CountType::Node, sequence.to_string()),
                name: gb.get_fname(),
                label: "Average Degree".to_string(),
                sequence: sequence.to_string(),
                values,
            })
            .collect();
        let table_text = self.generate_table(Some(gb))?;
        let table_text = format!("`{}`", table_text);
        let regional_degree_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{}", CountType::Node),
            analysis: "Regional Degree".to_string(),
            table: Some(table_text),
            run_name: self.get_run_name(gb),
            run_id: self.get_run_id(gb),
            countable: CountType::Node.to_string(),
            items,
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(regional_degree_tabs)
    }

    fn get_graph_requirements(&self) -> std::collections::HashSet<super::InputRequirement> {
        HashSet::from([
            InputRequirement::Node,
            InputRequirement::Bp,
            InputRequirement::Edge,
            InputRequirement::Hist,
            InputRequirement::Path(self.reference_text.clone()),
        ])
    }

    fn get_type(&self) -> String {
        "RegionalDegree".to_string()
    }
}

impl ConstructibleAnalysis for RegionalDegree {
    fn from_parameter(parameter: crate::analysis_parameter::AnalysisParameter) -> Self {
        let (reference, window_size) = match parameter {
            AnalysisParameter::RegionalDegree {
                reference,
                window_size,
            } => (reference, window_size),
            _ => panic!("Regional degree should only be called with correct parameter"),
        };
        Self {
            reference_text: PathSegment::from_str(&reference),
            window_size,
            values: HashMap::new(),
        }
    }
}

impl RegionalDegree {
    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-regionaldegree", gb.get_run_id())
    }

    fn set_values(&mut self, gb: &GraphBroker) {
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
        let degrees = gb.get_degree();
        let mut all_degrees_of_windows: HashMap<PathSegment, Vec<(f64, usize, usize)>> =
            HashMap::new();
        for (sequence_id, sequence) in ref_paths {
            for (contig_id, contig) in sequence {
                let windows = get_windows(contig, node_lens, self.window_size, &neighbors);
                let contig_start = contig_id.start.unwrap_or_default();
                let degrees_of_windows: Vec<(f64, usize, usize)> = windows
                    .iter()
                    .map(|(window, start, end)| {
                        (
                            window
                                .iter()
                                .map(|(node, _)| degrees[node.0 as usize] as f64)
                                .sum::<f64>()
                                / window.len() as f64,
                            *start + contig_start,
                            *end + contig_start,
                        )
                    })
                    .collect();
                all_degrees_of_windows
                    .entry(sequence_id.clone())
                    .or_default()
                    .extend(degrees_of_windows);
            }
        }
        self.values = all_degrees_of_windows;
    }
}
