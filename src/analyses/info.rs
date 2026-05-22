use itertools::Itertools;

use crate::{
    analyses::{AnalysisSection, MatrixBasedAnalysis},
    coverage_matrix::CoverageMatrix,
    html_report::ReportItem,
    util::get_default_plot_downloads,
};

pub struct Info {}

impl MatrixBasedAnalysis for Info {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        let mut res = format!(
            "# {}\n",
            std::env::args().collect::<Vec<String>>().join(" ")
        );
        res.push_str(matrix.get_file_info().to_string().as_str());
        Ok(res)
    }

    fn get_type(&self) -> String {
        "Info".to_string()
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);

        let info_values = matrix
            .get_file_info()
            .iterate_infos()
            .map(|(k, v)| vec![k.to_string(), v.to_string()])
            .collect_vec();

        let run_name = matrix.get_run_name();
        let run_id = matrix.get_run_id();
        let safe_run_name = run_id.to_lowercase().replace(&[' ', '|', '\\'], "-");
        Ok(vec![AnalysisSection {
            id: format!("{safe_run_name}-graph"),
            analysis: "Pangenome Info".to_string(),
            run_name: run_name.to_string(),
            run_id: run_id.to_string(),
            countable: "Graph Info".to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::Table {
                id: "info-1-table".to_string(),
                header: vec!["Statistic".to_string(), "Value".to_string()],
                values: info_values,
            }],
            plot_downloads: get_default_plot_downloads(),
        }])
    }
}

impl Info {
    pub fn new() -> Self {
        Self {}
    }
}

#[derive(Debug)]
pub struct FileInfo {
    infos: Vec<(String, String)>,
    filetype: String,
}

impl FileInfo {
    pub fn new(filetype: &str) -> Self {
        Self {
            infos: Vec::new(),
            filetype: filetype.to_string(),
        }
    }

    pub fn get_filetype(&self) -> &str {
        &self.filetype
    }

    pub fn add_info(&mut self, key: &str, value: &str) {
        self.infos.push((key.to_string(), value.to_string()));
    }

    pub fn iterate_infos(&self) -> impl Iterator<Item = (&str, &str)> + '_ {
        self.infos.iter().map(|(a, b)| (a.as_str(), b.as_str()))
    }

    #[cfg(test)]
    fn get_info_count(&self) -> usize {
        self.infos.len()
    }
}

impl std::fmt::Display for FileInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = self
            .infos
            .iter()
            .map(|(k, v)| {
                let k_c = cleanup_string_for_printing(k);
                let v_c = cleanup_string_for_printing(v);
                format!("{}\t{}", k_c, v_c)
            })
            .join("\n");
        write!(f, "{}", text)
    }
}

fn cleanup_string_for_printing(t: &str) -> String {
    t.to_lowercase().replace(|ch: char| ch.is_whitespace(), "-")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_simple_gfa_file_info() -> FileInfo {
        let mut fileinfo = FileInfo::new("gfa");
        fileinfo.add_info("Number of nodes", "12");
        fileinfo
    }

    #[test]
    fn test_file_info_creation() {
        let file_type = "gfa";
        let file_info = FileInfo::new(file_type);
        let returned_filetype = file_info.get_filetype();
        assert_eq!(file_type, returned_filetype);
    }

    #[test]
    fn test_file_info_addition() {
        let fileinfo = create_simple_gfa_file_info();
        assert_eq!(fileinfo.get_info_count(), 1);
    }

    #[test]
    fn test_file_info_printing() {
        let mut fileinfo = create_simple_gfa_file_info();
        fileinfo.add_info("Number of edges", "6");
        let text = "number-of-nodes\t12\nnumber-of-edges\t6";
        assert_eq!(&fileinfo.to_string(), text);
    }

    #[test]
    fn test_file_info_iterating() {
        let mut fileinfo = create_simple_gfa_file_info();
        fileinfo.add_info("Number of edges", "6");
        let mut iter = fileinfo.iterate_infos();
        assert_eq!(iter.next().unwrap(), ("Number of nodes", "12"));
        assert_eq!(iter.next().unwrap(), ("Number of edges", "6"));
    }
}
