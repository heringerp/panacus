use std::collections::HashSet;

use crate::{
    analyses::HistBasedAnalysis,
    hist::Hist,
    html_report::{AnalysisSection, ReportItem},
    io::write_table_with_start_index,
    util::{get_default_plot_downloads, CountType},
};

use super::InputRequirement;

pub struct CoverageLine {}

impl HistBasedAnalysis for CoverageLine {
    fn get_type(&self) -> String {
        "CoverageLine".to_string()
    }

    fn generate_table(&mut self, hist: &Hist) -> anyhow::Result<String> {
        log::info!("reporting coverage line table");
        let mut res = String::new();
        res.push_str(&crate::io::write_metadata_comments()?);

        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            String::new(),
            String::new(),
        ]];
        let mut output_columns = Vec::new();
        output_columns.push(
            hist.get_hist_values()
                .iter()
                .map(|x| *x as f64)
                .skip(1)
                .collect(),
        );
        header_cols.push(vec![
            "hist".to_string(),
            hist.get_feature_type().to_string(),
            String::new(),
            String::new(),
        ]);
        res.push_str(&write_table_with_start_index(
            &header_cols,
            &output_columns,
            1,
        )?);
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        hist: &Hist,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let table = self.generate_table(hist)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "coverage-line-{}",
            hist.get_run_id()
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let mut values: Vec<_> = hist.get_hist_values().clone();
        while let Some(last) = values.pop() {
            if last != 0 {
                values.push(0);
                break;
            }
        }
        let values: Vec<f32> = values.into_iter().skip(1).map(|c| c as f32).collect();
        let k = hist.get_feature_type();
        let coverage_line_tabs = AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Coverage Line".to_string(),
            table: Some(table.clone()),
            run_name: hist.get_run_name().to_owned(),
            run_id: hist.get_run_id().to_owned(),
            countable: k.to_string(),
            items: vec![ReportItem::Line {
                id: format!("{id_prefix}-{k}"),
                name: hist.get_run_name().to_owned(),
                x_label: "Allele count".to_string(),
                y_label: format!("#{}s", k),
                x_values: (1..=values.len()).map(|s| s as f32).collect(),
                y_values: values,
                log_x: true,
                log_y: true,
            }],
            plot_downloads: get_default_plot_downloads(),
        };
        Ok(vec![coverage_line_tabs])
    }
}

impl CoverageLine {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }

    pub fn new() -> Self {
        Self {}
    }
}
