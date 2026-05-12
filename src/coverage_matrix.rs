use crate::{file_formats::gfa_parser::SparseMatrix, hist::Hist, util::ItemTable};

pub struct CoverageMatrix {
    count_of_features: usize,
    feature_lengths: Vec<usize>,
    feature_positions: Vec<usize>,
    feature_names: Vec<String>,
    path_names: Vec<String>,
    matrix: SparseMatrix,
    feature_type: String,
    run_id: String,
    run_name: String,
}

impl CoverageMatrix {
    // TODO: Add FILE info object
    pub fn new(feature_type: String, run_id: String, run_name: String) -> Self {
        Self {
            count_of_features: 0,
            feature_lengths: Vec::new(),
            feature_positions: Vec::new(),
            feature_names: Vec::new(),
            path_names: Vec::new(),
            matrix: SparseMatrix::new(),
            feature_type,
            run_id,
            run_name,
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

    pub fn get_regional_iterator(&self, _window_size: usize, _window_step_size: usize) {
        unimplemented!()
    }

    pub fn get_feature_counts_for_paths(&self, _paths: &Vec<usize>) -> Vec<usize> {
        unimplemented!()
    }

    pub fn get_feature_counts_for_subset(
        &self,
        _features: &Vec<usize>,
        _paths: &Vec<usize>,
    ) -> Vec<usize> {
        unimplemented!()
    }

    pub fn get_feature_counts(&self) -> Vec<usize> {
        unimplemented!()
    }

    pub fn get_feature_names(&self) -> Vec<String> {
        self.feature_names.clone()
    }

    pub fn get_feature_lengths(&self) -> &Vec<usize> {
        &self.feature_lengths
    }

    pub fn get_counts_for_feature(&self, _id: usize) -> Vec<usize> {
        unimplemented!()
    }

    pub fn get_path_names(&self) -> &Vec<String> {
        &self.path_names
    }

    pub fn get_feature_count(&self) -> usize {
        self.count_of_features
    }

    pub fn get_feature_type(&self) -> &str {
        &self.feature_type
    }

    pub fn get_run_id(&self) -> &str {
        &self.run_id
    }

    pub fn get_run_name(&self) -> &str {
        &self.run_name
    }

    pub fn get_file_info(&self) -> &str {
        unimplemented!()
    }
}
