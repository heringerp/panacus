use std::ops::Index;

pub struct Hist {
    count_of_features: usize,
    hist: Vec<usize>,
    feature_type: String,
    run_id: String,
    run_name: String,
}

impl Hist {
    /// Creates a histogram with a given maximum coverage.
    pub fn from_maximum_coverage(
        maximum_coverage: usize,
        feature_type: String,
        run_id: String,
        run_name: String,
    ) -> Self {
        Self {
            count_of_features: 0,
            hist: vec![0; maximum_coverage + 1],
            feature_type,
            run_id,
            run_name,
        }
    }

    #[cfg(test)]
    pub fn from_values(values: Vec<usize>) -> Self {
        Self {
            count_of_features: 0,
            hist: values,
            feature_type: "Test".to_string(),
            run_id: "".to_string(),
            run_name: "".to_string(),
        }
    }

    /// Inserts a single feature with a given coverage into the histogram.
    /// Will panic if the coverage is bigger than the maximum coverage.
    pub fn insert_feature_of_coverage(&mut self, coverage: usize) {
        self.hist[coverage] += 1;
        self.count_of_features += 1;
    }

    pub fn insert_feature_of_coverage_and_length(&mut self, coverage: usize, length: usize) {
        self.hist[coverage] += length;
        self.count_of_features += length;
    }

    /// Gets the number of features that are currently inserted into
    /// the histogram (this is multiplied by the length of the feature)
    pub fn get_number_of_features(&self) -> usize {
        self.count_of_features
    }

    /// Gets the maximum allowed coverage for this histogram.
    pub fn get_maximum_coverage(&self) -> usize {
        self.hist.len() - 1
    }

    pub fn get_hist_values(&self) -> &[usize] {
        &self.hist
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
}

impl Index<usize> for Hist {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.hist[index]
    }
}
