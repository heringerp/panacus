use crate::{coverage_matrix::CoverageMatrix, hist::Hist};

pub mod gfa_parser;

pub trait FileFormatParser {
    fn generate_hist(self) -> Hist;
    fn generate_matrix(self) -> CoverageMatrix;
    fn get_file_info(&self);
}
