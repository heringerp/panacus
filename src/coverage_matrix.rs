use std::collections::HashMap;

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
        let mut hist = Hist::from_maximum_coverage(
            self.path_names.len(),
            self.feature_type.clone(),
            self.run_id.clone(),
            self.run_name.clone(),
        );
        self.get_feature_counts()
            .into_iter()
            .zip(self.feature_lengths.iter())
            .for_each(|(c, l)| hist.insert_feature_of_coverage_and_length(c, *l));
        hist
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
        feature_names: Vec<String>,
        item_table: ItemTable,
    ) {
        self.path_names = path_names;
        self.feature_lengths = feature_lengths;
        self.count_of_features = self.feature_lengths.len();
        self.feature_names = feature_names;
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

    pub fn get_feature_name(&self, feature: usize) -> &str {
        &self.feature_names[feature]
    }

    /// Creates an iterator over indices in the order (in order),
    /// containing only the indices for which there exists an element
    pub fn get_appearances_for_order(&self, feature: usize, order: &Vec<String>) -> Vec<usize> {
        let translation_table: HashMap<&String, usize> = self
            .path_names
            .iter()
            .enumerate()
            .map(|(idx, path_name)| (path_name, idx))
            .collect();
        order
            .iter()
            .enumerate()
            .filter_map(|(idx, path_name)| {
                let translated_idx = translation_table[path_name];
                if self.matrix.contains(feature, translated_idx as u64) {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_feature_counts(&self) -> Vec<usize> {
        (0..self.feature_lengths.len())
            .map(|i| self.matrix.get_feature_occurrence_count(i))
            .collect()
    }

    pub fn get_feature_names(&self) -> Vec<String> {
        self.feature_names.clone()
    }

    pub fn get_feature_lengths(&self) -> &Vec<usize> {
        &self.feature_lengths
    }

    pub fn get_count_of_feature(&self, feature: usize) -> usize {
        self.matrix.get_feature_occurrence_count(feature) as usize
    }

    pub fn get_counts_for_feature(&self, id: usize) -> Vec<usize> {
        self.matrix
            .get_counts_for_feature(id, self.path_names.len())
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
