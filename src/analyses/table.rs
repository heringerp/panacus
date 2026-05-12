use crate::{
    analyses::{InputRequirement, MatrixBasedAnalysis},
    coverage_matrix::CoverageMatrix,
    io::write_metadata_comments,
    util::CountType,
};
use std::collections::HashSet;

use super::AnalysisSection;

pub struct Table {
    total: bool,
    order: Option<String>,
}

impl MatrixBasedAnalysis for Table {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        let mut string = write_metadata_comments()?;
        string.push_str("panacus\ttable\n");
        string.push_str("feature");
        for path_name in matrix.get_path_names() {
            string.push_str("\t");
            string.push_str(&path_name);
        }
        string.push_str("\n");
        for i in 0..matrix.get_feature_count() {
            string.push_str(&i.to_string());
            for value in matrix.get_counts_for_feature(i) {
                string.push_str("\t");
                string.push_str(&format!("{}", value));
            }
            string.push_str("\n");
        }
        Ok(string)
    }

    fn get_type(&self) -> String {
        "Table".to_string()
    }

    fn generate_report_section(
        &mut self,
        _matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        Ok(Vec::new())
    }
}

impl Table {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }

    pub fn new(total: bool, order: Option<String>) -> Self {
        Self { total, order }
    }
}
