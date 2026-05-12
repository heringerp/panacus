use kodama::Dendrogram;

use crate::analyses::MatrixBasedAnalysis;
use crate::analysis_parameter::ClusterMethod;
use crate::coverage_matrix::CoverageMatrix;
use crate::util::get_default_plot_downloads;
use crate::{
    analyses::InputRequirement, html_report::ReportItem, io::write_metadata_comments,
    util::CountType,
};
use std::collections::HashSet;
use std::usize;

use super::AnalysisSection;

pub struct Similarity {
    cluster_method: ClusterMethod,
    table: Option<Vec<Vec<f32>>>,
    labels: Option<Vec<String>>,
}

impl MatrixBasedAnalysis for Similarity {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        if self.table.is_none() {
            self.set_table(matrix);
        }
        let mut text = write_metadata_comments()?;
        let table = self.table.as_ref().unwrap();
        let labels = self.labels.as_ref().unwrap();
        text.push_str(&get_table_string(table, labels));
        Ok(text)
    }

    fn get_type(&self) -> String {
        "Similarity".to_string()
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        if self.table.is_none() {
            self.set_table(matrix);
        }

        let k = matrix.get_feature_type();
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "sim-heat-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Similarity Heatmap".to_string(),
            table: Some(table.clone()),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: k.to_string(),
            items: vec![ReportItem::Heatmap {
                id: format!("{id_prefix}-{k}"),
                name: matrix.get_run_name().to_owned(),
                x_labels: self.labels.as_ref().unwrap().clone(),
                y_labels: self.labels.as_ref().unwrap().clone(),
                values: self.table.as_ref().unwrap().clone(),
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(tabs)
    }
}

impl Similarity {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }

    pub fn new(cluster_method: ClusterMethod) -> Self {
        Self {
            cluster_method,
            table: None,
            labels: None,
        }
    }

    fn set_table(&mut self, _matrix: &CoverageMatrix) {
        unimplemented!()
    }
}

fn sort_by_indices<T>(list: &mut Vec<T>, indices: &Vec<usize>) {
    let mut indices = indices.clone();
    for i in 0..indices.len() {
        while i != indices[i] {
            let new_i = indices[i];
            indices.swap(i, new_i);
            list.swap(i, new_i);
        }
    }
}

fn get_order_from_dendrogram(dend: &Dendrogram<f32>) -> Vec<usize> {
    let observations = dend.observations();
    let mut indices = Vec::new();
    for step in dend.steps() {
        if step.cluster1 < observations {
            indices.push(step.cluster1);
        }
        if step.cluster2 < observations {
            indices.push(step.cluster2);
        }
    }
    indices
}

fn get_table_string(table: &Vec<Vec<f32>>, groups: &Vec<String>) -> String {
    let mut res = String::new();
    res.push_str("group");
    for group in groups {
        res.push_str(&format!("\t{}", group));
    }
    res.push_str("\n");
    for (row_index, row) in table.iter().enumerate() {
        res.push_str(&groups[row_index]);
        for cell in row {
            res.push_str(&format!("\t{}", cell));
        }
        res.push_str("\n");
    }
    res
}

fn euclidean(row1: &Vec<f32>, row2: &Vec<f32>) -> f32 {
    row1.iter()
        .zip(row2.iter())
        .map(|(v1, v2)| (v1 - v2).powf(2.0))
        .sum::<f32>()
        .sqrt()
}

fn calculate_distances(table: &Vec<Vec<f32>>) -> Vec<f32> {
    let mut condensed = vec![];
    for row in 0..table.len() - 1 {
        for col in row + 1..table.len() {
            condensed.push(euclidean(&table[row], &table[col]));
        }
    }
    condensed
}
