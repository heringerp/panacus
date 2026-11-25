use crate::analyses::utils::growth::calc_all_growths;
use crate::analyses::HistBasedAnalysis;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::hist::Hist;
use crate::html_report::ReportItem;
use crate::io::write_table;
use crate::util::get_default_plot_downloads;

use super::AnalysisSection;

pub struct Growth {
    coverage: Option<String>,
    quorum: Option<String>,
    add_hist: bool,
    growths: Option<Vec<Vec<f64>>>,
    thresholds: Option<ThresholdContainer>,
}

impl HistBasedAnalysis for Growth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(&mut self, hist: &Hist) -> anyhow::Result<String> {
        log::info!("reporting hist table, {}", self.add_hist,);

        self.set_inner(hist)?;
        let growths = self.growths.as_ref().unwrap();
        let thresholds = self.thresholds.as_ref().unwrap();

        let mut res = String::new();
        res.push_str(&format!(
            "# {}\n",
            std::env::args().collect::<Vec<String>>().join(" ")
        ));
        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            "coverage".to_string(),
            "quorum".to_string(),
        ]];
        let mut output_columns: Vec<Vec<f64>> = Vec::new();

        let hist_values = hist.get_hist_values();

        if self.add_hist {
            output_columns.push(hist_values.iter().map(|x| *x as f64).collect());
            header_cols.push(vec![
                "hist".to_string(),
                hist.get_feature_type().to_string(),
                String::new(),
                String::new(),
            ])
        }

        let growths: Vec<Vec<f64>> = growths
            .iter()
            .map(|row| {
                row.iter()
                    .map(|el| if el.is_nan() { 0.0 } else { *el })
                    .collect()
            })
            .collect();
        output_columns.extend(growths);
        let m = thresholds.coverage.len();
        header_cols.extend(
            std::iter::repeat_n("growth", m)
                .zip(std::iter::repeat_n(hist.get_feature_type(), m))
                .zip(&thresholds.coverage)
                .zip(&thresholds.quorum)
                .map(|(((p, t), c), q)| {
                    vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                }),
        );
        res.push_str(&write_table(&header_cols, &output_columns)?);
        Ok(res)
    }

    fn generate_report_section(&mut self, hist: &Hist) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(hist)?;
        let hist_aux = &self.thresholds.as_ref().unwrap();
        let growth_labels = (0..hist_aux.coverage.len())
            .map(|i| {
                format!(
                    "coverage ≥ {}, quorum ≥ {}%",
                    hist_aux.coverage[i].get_string(),
                    match hist_aux.quorum[i] {
                        crate::util::Threshold::Relative(x) => (x * 100.0).to_string(),
                        crate::util::Threshold::Absolute(x) => (x * 100).to_string(),
                    }
                )
            })
            .collect::<Vec<_>>();
        let table = self.generate_table(hist)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "pan-growth-{}",
            hist.get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        let k = hist.get_feature_type();
        let growths = self.growths.as_ref().unwrap();
        let growth_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Pangenome Growth".to_string(),
            run_name: hist.get_run_name().to_owned(),
            run_id: hist.get_run_id().to_owned(),
            countable: hist.get_feature_type().to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::MultiBar {
                id: format!("{id_prefix}-{k}"),
                names: growth_labels.clone(),
                x_label: "taxa".to_string(),
                y_label: format!("#{}s", k),
                labels: (0..(growths[0].len())).map(|i| i.to_string()).collect(),
                values: growths
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|el| if el.is_nan() { 0.0 } else { *el })
                            .collect()
                    })
                    .collect(),
                curve: None,
                alpha: None,
                log_toggle: false,
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(growth_tabs)
    }
}

impl Growth {
    pub fn new(coverage: Option<String>, quorum: Option<String>, add_hist: bool) -> Self {
        Self {
            coverage,
            quorum,
            add_hist,
            growths: None,
            thresholds: None,
        }
    }

    fn set_inner(&mut self, hist: &Hist) -> anyhow::Result<()> {
        if self.growths.is_some() {
            return Ok(());
        }
        let quorum = self.quorum.to_owned().unwrap_or("0".to_string());
        let coverage = self.coverage.to_owned().unwrap_or("1".to_string());
        let hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;

        let growths = calc_all_growths(hist, &hist_aux, true);
        self.growths = Some(growths);
        self.thresholds = Some(hist_aux);
        Ok(())
    }
}
