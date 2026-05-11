use crate::{file_formats::gfa_parser::SparseMatrix, hist::Hist, util::ItemTable};

pub struct CoverageMatrix {
    count_of_features: usize,
    feature_lengths: Vec<usize>,
    feature_positions: Vec<usize>,
    path_names: Vec<String>,
    matrix: SparseMatrix,
}

impl CoverageMatrix {
    pub fn new() -> Self {
        Self {
            count_of_features: 0,
            feature_lengths: Vec::new(),
            feature_positions: Vec::new(),
            path_names: Vec::new(),
            matrix: SparseMatrix::new(),
        }
    }

    /// Calculates the histogram from the matrix. This should
    /// only be used if the CoverageMatrix is needed anyways.
    /// If it isn't needed generate the hist directly from the
    /// file (using `FileFormatParser`'s `generate_hist()`).
    pub fn get_hist(&self) -> Hist {
        unimplemented!()
    }

    pub fn insert_feature(
        &mut self,
        feature_length: usize,
        feature_position: usize,
        feature: Vec<u32>,
    ) {
        self.count_of_features += 1;
        self.feature_lengths.push(feature_length);
        self.feature_positions.push(feature_position);
        self.matrix.insert_row(feature);
    }

    pub fn insert_item_table(
        &mut self,
        path_names: Vec<String>,
        feature_lengths: Vec<usize>,
        feature_positions: Vec<usize>,
        item_table: ItemTable,
    ) {
        self.path_names = path_names;
        self.count_of_features = feature_lengths.len();
        self.feature_lengths = feature_lengths;
        self.feature_positions = feature_positions;
        self.matrix
            .insert_item_table(self.count_of_features, item_table);
    }

    pub fn get_regional_iterator(&self, window_size: usize, window_step_size: usize) {
        unimplemented!()
    }

    pub fn get_feature_counts_for_paths(&self, paths: &Vec<usize>) -> Vec<usize> {
        unimplemented!()
    }

    pub fn get_feature_counts_for_subset(
        &self,
        features: &Vec<usize>,
        paths: &Vec<usize>,
    ) -> Vec<usize> {
        unimplemented!()
    }
}
