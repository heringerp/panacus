use std::collections::{HashMap, HashSet};

use crate::{
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, PathSegment},
    html_report::{AnalysisSection, ReportItem},
    util::get_default_plot_downloads,
};

use super::{
    regional_helpers::{get_close_nodes, get_ref_length, get_windows},
    Analysis, ConstructibleAnalysis, InputRequirement,
};

const NUMBER_OF_WINDOWS: usize = 1000;

pub struct RegionalDegree {
    reference: PathSegment,
    window_size: usize,
    values: Vec<f64>,
    ref_len: usize,
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
        for (i, entry) in self.values.iter().enumerate() {
            let line = format!(
                "{} {} {} {}\n",
                self.reference.to_string(),
                i * self.window_size,
                if i >= self.values.len() - 1 {
                    self.ref_len - 1
                } else {
                    (i + 1) * self.window_size - 1
                },
                *entry
            );
            text.push_str(&line);
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
        let values: Vec<(usize, usize, f64)> = self
            .values
            .iter()
            .enumerate()
            .map(|(i, v)| {
                if i < self.values.len() - 2 {
                    (i * self.window_size, (i + 2) * self.window_size, *v)
                } else {
                    (i * self.window_size, self.ref_len - 1, *v)
                }
            })
            .collect();
        let values = HashMap::from([(self.reference.to_string(), values)]);
        let id_prefix = format!(
            "regional-degree-{}",
            self.get_run_id(gb)
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let table_text = self.generate_table(Some(gb))?;
        let table_text = format!("`{}`", table_text);
        let regional_degree_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}"),
            analysis: "Regional".to_string(),
            table: Some(table_text),
            run_name: self.get_run_name(gb),
            run_id: self.get_run_id(gb),
            countable: "Degree".to_string(),
            items: vec![ReportItem::Chromosomal {
                id: format!("{id_prefix}"),
                name: gb.get_fname(),
                label: "Average Degree".to_string(),
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
            window_size: 0,
            values: Vec::new(),
            ref_len: 0,
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

    fn set_values(&mut self, gb: &GraphBroker) {
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
        self.window_size = window_size;
        self.values = degrees_of_windows;
        self.ref_len = get_ref_length(ref_nodes, node_lens) as usize;
    }
}
