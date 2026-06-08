pub mod coverage_colors;
pub mod coverage_line;
pub mod growth;
pub mod hist;
pub mod info;
pub mod node_distribution;
pub mod ordered_histgrowth;
pub mod regional_growth;
pub mod regional_helpers;
pub mod regional_variation;
pub mod section_growth;
pub mod similarity;
pub mod table;

use std::collections::HashSet;

use crate::{
    coverage_matrix::CoverageMatrix, file_formats::gfa_parser::PathSegment, hist::Hist,
    html_report::AnalysisSection, util::CountType,
};

pub trait HistBasedAnalysis {
    fn generate_table(&mut self, hist: &Hist) -> anyhow::Result<String>;
    fn generate_report_section(&mut self, hist: &Hist) -> anyhow::Result<Vec<AnalysisSection>>;
    fn get_type(&self) -> String;
}

pub trait MatrixBasedAnalysis {
    fn generate_table(&mut self, hist: &CoverageMatrix) -> anyhow::Result<String>;
    fn generate_report_section(
        &mut self,
        hist: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>>;
    fn get_type(&self) -> String;
}

pub enum Analysis {
    HistBased(Box<dyn HistBasedAnalysis>),
    MatrixBased(Box<dyn MatrixBasedAnalysis>),
}

#[derive(PartialEq, Eq, Debug, Clone, Hash)]
pub enum InputRequirement {
    Node,
    Edge,
    Bp,
    PathLens,
    Hist,
    AbacusByGroup(CountType),
    Graph(String),
    Path(PathSegment),
}

impl InputRequirement {
    pub fn from_count(count: CountType) -> HashSet<Self> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }
}
