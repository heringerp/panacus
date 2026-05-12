use std::collections::{HashMap, HashSet};

use crate::{
    analyses::{AnalysisSection, MatrixBasedAnalysis},
    coverage_matrix::CoverageMatrix,
    file_formats::gfa_parser::{Edge, ItemId},
};

pub struct Info {}

impl MatrixBasedAnalysis for Info {
    fn generate_table(&mut self, _matrix: &CoverageMatrix) -> anyhow::Result<String> {
        unimplemented!()
    }

    fn get_type(&self) -> String {
        "Info".to_string()
    }

    fn generate_report_section(
        &mut self,
        _matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        unimplemented!()
    }
}

impl Info {
    fn set_info(&mut self, _matrix: &CoverageMatrix) {
        unimplemented!()
    }
}

fn connected_components(edge2id: &HashMap<Edge, ItemId>, nodes: &Vec<ItemId>) -> Vec<u32> {
    let mut component_lengths = Vec::new();
    let mut visited: HashSet<ItemId> = HashSet::new();
    let edges: HashMap<ItemId, Vec<ItemId>> = edge2id
        .keys()
        .map(|x| (x.0, x.2))
        .chain(edge2id.keys().map(|x| (x.2, x.0)))
        .fold(HashMap::new(), |mut acc, (k, v)| {
            acc.entry(k).and_modify(|x| x.push(v)).or_insert(vec![v]);
            acc
        });
    for node in nodes {
        if !visited.contains(node) {
            component_lengths.push(dfs(&edges, *node, &mut visited));
        }
    }
    component_lengths
}

fn dfs(edges: &HashMap<ItemId, Vec<ItemId>>, node: ItemId, visited: &mut HashSet<ItemId>) -> u32 {
    let mut s = Vec::new();
    let mut length = 0;
    s.push(node);
    while let Some(v) = s.pop() {
        if visited.contains(&v) {
            continue;
        }
        visited.insert(v);
        length += 1;
        if !edges.contains_key(&v) {
            continue;
        }
        for neigh in &edges[&v] {
            if !visited.contains(neigh) {
                s.push(*neigh);
            }
        }
    }
    length
}
