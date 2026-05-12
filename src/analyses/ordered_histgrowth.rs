use std::collections::HashSet;

use crate::analyses::InputRequirement;
use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::html_report::ReportItem;
use crate::util::{get_default_plot_downloads, CountType};

use super::AnalysisSection;

type Growths = Vec<Vec<f64>>;

pub struct OrderedHistgrowth {
    coverage: Option<String>,
    quorum: Option<String>,
    order: Option<String>,
    inner: Option<InnerOrderedGrowth>,
}

impl MatrixBasedAnalysis for OrderedHistgrowth {
    fn get_type(&self) -> String {
        "OrderedHistgrowth".to_string()
    }
    fn generate_table(&mut self, _matrix: &CoverageMatrix) -> anyhow::Result<String> {
        unimplemented!()
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(matrix)?;
        let hist_aux = &self.inner.as_ref().unwrap().hist_aux;
        let growth_labels = (0..hist_aux.coverage.len())
            .map(|i| {
                format!(
                    "coverage ≥ {}, quorum ≥ {}%",
                    hist_aux.coverage[i].get_string(),
                    hist_aux.quorum[i].get_string()
                )
            })
            .collect::<Vec<_>>();
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let growths = &self.inner.as_ref().unwrap().growths;
        let id_prefix = format!(
            "pan-ordered-growth-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let labels = matrix.get_path_names().clone();
        let count = matrix.get_feature_type();
        let growth_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}"),
            analysis: "Ordered Growth".to_string(),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: matrix.get_feature_type().to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::MultiBar {
                id: format!("{id_prefix}"),
                names: growth_labels.clone(),
                x_label: "taxa".to_string(),
                y_label: format!("{}s", count),
                //labels: (1..growths[0].len()).map(|i| i.to_string()).collect(),
                labels,
                values: growths.clone(),
                log_toggle: false,
                alpha: None,
                curve: None,
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(growth_tabs)
    }
}

impl OrderedHistgrowth {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }

    fn set_inner(&mut self, _matrix: &CoverageMatrix) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        let quorum = self.quorum.clone().unwrap_or("0".to_string());
        let coverage = self.coverage.clone().unwrap_or("1".to_string());
        let _hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;

        unimplemented!();
    }
}

struct InnerOrderedGrowth {
    growths: Growths,
    hist_aux: ThresholdContainer,
}
