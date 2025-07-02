use std::{
    cell::Ref,
    collections::{HashMap, HashSet, VecDeque},
};

use crate::{
    analysis_parameter::AnalysisParameter,
    graph_broker::{GraphBroker, ItemId, PathSegment},
    html_report::{AnalysisSection, ReportItem},
    util::get_default_plot_downloads,
};

use super::{Analysis, ConstructibleAnalysis, InputRequirement};

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
        let close_ones_of_ref = Self::get_close_nodes(ref_nodes, &neighbors);
        // eprintln!("paths: {:?}", gb.get_path(&self.reference));
        let node_lens = gb.get_node_lens();
        let (windows, window_size) = Self::get_windows(ref_nodes, node_lens);
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
    fn get_close_nodes(
        ref_nodes: &Vec<ItemId>,
        neighbors: &HashMap<ItemId, Vec<ItemId>>,
    ) -> HashMap<ItemId, Vec<ItemId>> {
        let mut queue: VecDeque<ItemId> = VecDeque::new();
        let mut closest_ref: HashMap<ItemId, ItemId> = HashMap::new();
        let ref_node_set: HashSet<ItemId> = ref_nodes.iter().copied().collect();
        for ref_node in ref_nodes {
            closest_ref.insert(*ref_node, *ref_node);
            for neighbor in neighbors.get(ref_node).unwrap_or(&Vec::new()) {
                if !ref_node_set.contains(neighbor) && !closest_ref.contains_key(neighbor) {
                    closest_ref.insert(*neighbor, *ref_node);
                    queue.push_front(*neighbor);
                }
            }
        }
        while let Some(current_node) = queue.pop_back() {
            for neighbor in neighbors.get(&current_node).unwrap_or(&Vec::new()) {
                if !closest_ref.contains_key(neighbor) {
                    closest_ref.insert(*neighbor, current_node);
                    queue.push_front(*neighbor);
                }
            }
        }
        let mut flipped: HashMap<ItemId, Vec<ItemId>> = HashMap::new();
        for (entry_node, ref_node) in closest_ref {
            flipped.entry(ref_node).or_default().push(entry_node);
        }
        flipped
    }

    fn get_ref_length(ref_nodes: &Vec<ItemId>, node_lens: &Vec<u32>) -> u32 {
        let length = ref_nodes
            .iter()
            .map(|node| node_lens[node.0 as usize])
            .sum();
        length
    }

    fn get_windows(
        ref_nodes: &Vec<ItemId>,
        node_lens: &Vec<u32>,
    ) -> (Vec<Vec<(ItemId, usize)>>, usize) {
        let ref_length = Self::get_ref_length(ref_nodes, node_lens) as usize;
        let window_size = ref_length as usize / NUMBER_OF_WINDOWS;
        let window_size = if window_size * NUMBER_OF_WINDOWS == ref_length as usize {
            window_size
        } else {
            window_size + 1 // Add 1 to make last window uncomplete/total number stay at NUMBER_OF_WINDOWS
        };
        let mut windows: Vec<Vec<(ItemId, usize)>> = vec![Vec::new(); NUMBER_OF_WINDOWS];
        let mut bp_counter = 0;
        let mut current_window_index = 0;
        let mut current_node_index = 0;
        let mut already_used_bps_of_node = 0;
        while bp_counter < ref_length {
            let current_window_length: usize =
                windows[current_window_index].iter().map(|(_, l)| l).sum();
            let remaining_bps_in_window = window_size - current_window_length as usize;
            let current_node_length = node_lens[ref_nodes[current_node_index].0 as usize] as usize
                - already_used_bps_of_node;
            if current_node_length <= remaining_bps_in_window {
                windows[current_window_index]
                    .push((ref_nodes[current_node_index], current_node_length as usize));
                bp_counter += current_node_length;
                current_node_index += 1;
                already_used_bps_of_node = 0;
                if current_node_length == remaining_bps_in_window {
                    current_window_index += 1;
                }
            } else {
                let cut_node_length = remaining_bps_in_window;
                windows[current_window_index]
                    .push((ref_nodes[current_node_index], cut_node_length as usize));
                bp_counter += cut_node_length;
                current_window_index += 1;
                already_used_bps_of_node += cut_node_length;
            };
        }
        (windows, window_size)
    }

    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-coverageline", gb.get_run_id())
    }
}
