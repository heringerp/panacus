use std::collections::HashMap;

use crate::{
    analyses::info::FileInfo,
    file_formats::gfa_parser::SparseMatrix,
    hist::Hist,
    util::{ItemTable, Threshold},
};

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
    file_info: FileInfo,
}

impl CoverageMatrix {
    // TODO: Add FILE info object
    pub fn new(
        feature_type: String,
        run_id: String,
        run_name: String,
        file_info: FileInfo,
    ) -> Self {
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
            file_info,
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

    pub fn get_hist_for_reference(&self, reference: &str) -> Hist {
        let mut hist = Hist::from_maximum_coverage(
            self.path_names.len() - 1,
            self.feature_type.clone(),
            self.run_id.clone(),
            self.run_name.clone(),
        );
        let r_idx = self
            .path_names
            .iter()
            .enumerate()
            .filter_map(|(i, p)| if p == reference { Some(i) } else { None })
            .next()
            .expect("Reference is part of paths");
        self.get_feature_counts()
            .into_iter()
            .enumerate()
            .zip(self.feature_lengths.iter())
            .for_each(|((i, c), l)| {
                // Insert the feature only if it is not part of the reference
                if !self.matrix.contains(i, r_idx as u64) {
                    hist.insert_feature_of_coverage_and_length(c, *l)
                }
            });
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

    /// Returns a histogram and the features
    /// that really appear in the histogram (i.e. are non-zero).
    /// If features is Some(list), the list is used to subset the features
    pub fn get_hist_and_used_features(
        &self,
        features: Option<&Vec<usize>>,
        paths: &Vec<usize>,
    ) -> (Hist, Vec<usize>) {
        let c = Threshold::Absolute(1);
        self.get_hist_and_used_features_with_coverage(features, paths, &c)
    }

    pub fn get_hist_and_used_features_with_coverage(
        &self,
        features: Option<&Vec<usize>>,
        paths: &Vec<usize>,
        c: &Threshold,
    ) -> (Hist, Vec<usize>) {
        let (abacus, used_features) =
            self.get_abacus_and_used_features_with_coverage(features, paths, c);
        log::info!("Hist-needed abacus: {:?}", abacus);
        let mut hist = Hist::from_maximum_coverage(
            paths.len(),
            self.feature_type.to_string(),
            self.run_id.to_string(),
            self.run_name.to_string(),
        );
        for (i, c) in abacus.iter().enumerate() {
            hist.insert_feature_of_coverage_and_length(*c as usize, self.feature_lengths[i]);
        }
        (hist, used_features)
    }

    pub fn get_abacus_and_used_features_with_coverage(
        &self,
        features: Option<&Vec<usize>>,
        paths: &Vec<usize>,
        c: &Threshold,
    ) -> (Vec<usize>, Vec<usize>) {
        let all_features = (0..self.count_of_features).collect();
        let features = match features {
            Some(f) => f,
            None => &all_features,
        };
        let mut paths = paths.clone();
        paths.sort();
        let mut non_zeroes = Vec::new();
        let mut abacus = vec![0; self.count_of_features];
        for feature in features {
            let mut path_idx = 0;
            let mut count = 0;
            let occurrences = self.matrix.get_occurrences(*feature);
            // log::info!(
            //     "Looking at feature: {} | {:?}, {:?}",
            //     *feature,
            //     occurrences,
            //     paths
            // );
            for o in occurrences {
                let o = *o as usize;
                while path_idx < paths.len() && paths[path_idx] < o {
                    path_idx += 1;
                }
                if path_idx >= paths.len() {
                    break;
                } else if o < paths[path_idx] {
                    continue;
                } else if o == paths[path_idx] {
                    count += 1;
                    path_idx += 1;
                }
            }
            // log::info!("Count: {}", count);
            if count >= c.to_absolute(self.count_of_features) {
                non_zeroes.push(*feature);
                abacus[*feature] = count;
            }
        }
        (abacus, non_zeroes)
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

    pub fn get_features_for_paths(&self, _paths: &Vec<usize>) -> Vec<usize> {
        unimplemented!()
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

    pub fn get_counts_for_feature_in_order(&self, id: usize, order: &Vec<String>) -> Vec<usize> {
        let translation_table: HashMap<&String, usize> = self
            .path_names
            .iter()
            .enumerate()
            .map(|(idx, path_name)| (path_name, idx))
            .collect();
        order
            .iter()
            .map(|path_name| {
                let translated_idx = translation_table[path_name];
                if self.matrix.contains(id, translated_idx as u64) {
                    1
                } else {
                    0
                }
            })
            .collect()
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

    pub fn get_file_info(&self) -> &FileInfo {
        &self.file_info
    }
}
