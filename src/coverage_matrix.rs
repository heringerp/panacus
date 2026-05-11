use crate::{graph_broker::SparseMatrix, hist::Hist};

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
        feature: Vec<usize>,
    ) {
        self.count_of_features += 1;
        self.feature_lengths.push(feature_length);
        self.feature_positions.push(feature_position);
        // TODO Insert feature itself
        unimplemented!()
    }

    pub fn insert_path(&mut self, path_name: &str, path: Vec<usize>) {
        self.path_names.push(path_name.to_owned());
        // TODO Insert path itself
        unimplemented!()
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
