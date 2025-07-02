use std::collections::{HashMap, HashSet};

use crate::{
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, PathSegment},
    html_report::{AnalysisSection, ReportItem},
    util::get_default_plot_downloads,
};

use super::{
    regional_helpers::{get_close_nodes, get_windows},
    Analysis, ConstructibleAnalysis, InputRequirement,
};

const NUMBER_OF_WINDOWS: usize = 1000;

pub struct RegionalDegree {
    reference: PathSegment,
}

impl Analysis for RegionalDegree {
    fn generate_table(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        Ok("".to_string())
    }

    fn generate_report_section(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let gb = gb.expect("Regional degree should always have a graph");
        let edge2id = gb.get_edges();
        let neighbors: HashMap<ItemId, Vec<ItemId>> = edge2id
            .keys()
            .map(|x| (x.0, x.2))
            .chain(edge2id.keys().map(|x| (x.2, x.0)))
            .fold(HashMap::new(), |mut acc, (k, v)| {
                acc.entry(k).and_modify(|x| x.push(v)).or_insert(vec![v]);
                acc
            });
        let ref_nodes = gb.get_path(&self.reference);
        let close_ones_of_ref = get_close_nodes(ref_nodes, &neighbors);
        let node_lens = gb.get_node_lens();
        let (windows, window_size) = get_windows(ref_nodes, node_lens, NUMBER_OF_WINDOWS);
        let degrees = gb.get_degree();
        let degrees_of_windows: Vec<f64> = windows
            .iter()
            .map(|window| {
                window
                    .iter()
                    .map(|(node, l)| {
                        *l as f64
                            * (close_ones_of_ref[node]
                                .iter()
                                .map(|current_node| degrees[current_node.0 as usize])
                                .sum::<u32>() as f64
                                / close_ones_of_ref[node].len() as f64)
                    })
                    .sum::<f64>()
                    / window.iter().map(|(_, l)| l).sum::<usize>() as f64
            })
            .collect();
        let values: Vec<(usize, f64)> = degrees_of_windows
            .into_iter()
            .enumerate()
            .map(|(i, v)| (i * window_size, v))
            .collect();

        let id_prefix = format!(
            "regional-degree-{}",
            self.get_run_id(gb)
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let regional_degree_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}"),
            analysis: "Regional".to_string(),
            table: None,
            run_name: self.get_run_name(gb),
            run_id: self.get_run_id(gb),
            countable: "Degree".to_string(),
            items: vec![ReportItem::Chromosomal {
                id: format!("{id_prefix}"),
                name: gb.get_fname(),
                values,
            }],
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
            InputRequirement::Path(self.reference.clone()),
        ])
    }

    fn get_type(&self) -> String {
        "RegionalDegree".to_string()
    }
}

impl ConstructibleAnalysis for RegionalDegree {
    fn from_parameter(parameter: crate::analysis_parameter::AnalysisParameter) -> Self {
        let reference = match parameter {
            AnalysisParameter::RegionalDegree { reference } => reference,
            _ => panic!("Regional degree should only be called with correct parameter"),
        };
        Self {
            reference: PathSegment::from_str(&reference),
        }
    }
}

impl RegionalDegree {
    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-coverageline", gb.get_run_id())
    }
}
