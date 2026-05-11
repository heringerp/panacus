use crate::{
    file_formats::FileFormatParser,
    graph_broker::{Hist, SparseMatrix},
};

pub struct GfaParser {
    filename: String,
}

impl GfaParser {
    pub fn new(filename: &str) -> Self {
        Self {
            filename: filename.to_owned(),
        }
    }
}

impl FileFormatParser for GfaParser {
    fn generate_hist(self) -> Hist {
        unimplemented!()
    }

    fn generate_matrix(self) -> SparseMatrix {
        unimplemented!()
    }

    fn get_file_info(&self) {
        unimplemented!()
    }
}
