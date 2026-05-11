use crate::graph_broker::{Hist, SparseMatrix};

pub mod gfa_parser;

pub trait FileFormatParser {
    fn generate_hist(self) -> Hist;
    fn generate_matrix(self) -> SparseMatrix;
    fn get_file_info(&self);
}
