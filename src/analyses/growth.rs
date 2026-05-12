use core::str;
use std::cmp;

use ml_helpers::linear_regression::huber_regressor::{solve, HuberRegressor};

use crate::analyses::HistBasedAnalysis;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::hist::Hist;
use crate::html_report::ReportItem;
use crate::util::get_default_plot_downloads;

use super::AnalysisSection;

pub struct Growth {
    coverage: Option<String>,
    quorum: Option<String>,
    add_hist: bool,
    add_alpha: bool,
    inner: Option<InnerGrowth>,
}

impl HistBasedAnalysis for Growth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(&mut self, hist: &Hist) -> anyhow::Result<String> {
        log::info!(
            "reporting hist table, {}, {}",
            self.add_hist,
            self.add_alpha
        );

        self.set_inner(hist)?;
        let mut _res = String::new();

        unimplemented!()
    }

    fn generate_report_section(&mut self, hist: &Hist) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(hist)?;
        // let hist_aux = &self.inner.as_ref().unwrap().hist_aux;
        // let growth_labels = (0..hist_aux.coverage.len())
        //     .map(|i| {
        //         format!(
        //             "coverage ≥ {}, quorum ≥ {}%",
        //             hist_aux.coverage[i].get_string(),
        //             match hist_aux.quorum[i] {
        //                 crate::util::Threshold::Relative(x) => (x * 100.0).to_string(),
        //                 crate::util::Threshold::Absolute(x) => (x * 100).to_string(),
        //             }
        //         )
        //     })
        //     .collect::<Vec<_>>();
        let table = self.generate_table(hist)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "pan-growth-{}",
            hist.get_run_id()
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let k = hist.get_feature_type();
        let growth_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Pangenome Growth".to_string(),
            run_name: hist.get_run_name().to_owned(),
            run_id: hist.get_run_id().to_owned(),
            countable: hist.get_feature_type().to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::MultiBar {
                id: format!("{id_prefix}-{k}"),
                names: Vec::new(),
                x_label: "taxa".to_string(),
                y_label: format!("#{}s", k),
                labels: Vec::new(),
                values: Vec::new(),
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
    pub fn new(
        coverage: Option<String>,
        quorum: Option<String>,
        add_hist: bool,
        add_alpha: bool,
    ) -> Self {
        Self {
            coverage,
            quorum,
            add_hist,
            add_alpha,
            inner: None,
        }
    }

    pub fn generate_table_from_hist(&mut self, _file: &str) -> anyhow::Result<String> {
        unimplemented!()
    }

    fn set_inner(&mut self, _hist: &Hist) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        let quorum = self.quorum.to_owned().unwrap_or("0".to_string());
        let coverage = self.coverage.to_owned().unwrap_or("1".to_string());
        let _hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;

        unimplemented!();
        // let hists = hist.get_hist_values();
        // let heaps_curves = hist_aux.has_full_growth_at_idx().map(|index| {
        //     log::info!("Calculating heaps law");
        //     let heaps_curves: Vec<_> = growths
        //         .iter()
        //         .zip(hists.iter())
        //         .map(|((_count_type, growth), (_count_type2, hist))| {
        //             let growth = growth[index].clone();
        //             let growth_len = growth.len();
        //             let growth_last = *growth.last().unwrap();
        //             let hist: Vec<f64> = hist.coverage.iter().map(|x| *x as f64).collect();
        //             let x1 = 2.0f64;
        //             let y1 = growth[1];
        //             let x2 = growth.len() as f64 - 1.0;
        //             let y2 = growth_last;
        //             let (alpha, _offset) = get_regression(&hist);
        //             if alpha >= 10.0 {
        //                 // TODO change 10 back to 1
        //                 (alpha, None)
        //             } else {
        //                 let gamma = 1.0 - alpha;
        //                 let k = (y1 - y2) / (x1.powf(gamma) - x2.powf(gamma));
        //                 let c = y1 - k * x1.powf(gamma);
        //                 let curve_values = (1..=growth_len)
        //                     .map(|x| (x as f64).powf(gamma) * k + c)
        //                     .collect::<Vec<_>>();
        //                 (alpha, Some(curve_values))
        //             }
        //         })
        //         .collect();
        //     heaps_curves
        // });
    }
}

#[allow(dead_code)]
fn get_regression(hist: &Vec<f64>) -> (f64, f64) {
    let x: Vec<f64> = (1..hist.len()).map(|x| x as f64).collect();
    let log_x: Vec<f64> = x
        .iter()
        .map(|x| x.log10())
        .map(|x| {
            if x.is_infinite() && x.is_sign_negative() {
                -100000.0
            } else {
                x
            }
        })
        .collect();
    let y = hist.iter().skip(1).copied().collect::<Vec<f64>>();
    let log_y: Vec<f64> = y
        .iter()
        .map(|y| y.log10())
        .map(|y| {
            if y.is_infinite() && y.is_sign_negative() {
                -100000.0
            } else {
                y
            }
        })
        .collect();
    let n = cmp::max(hist.len() / 2, 10);
    let log_x2: Vec<f64> = log_x.iter().skip(1).take(n).copied().collect();
    let log_y2: Vec<f64> = log_y.iter().skip(1).take(n).copied().collect();
    let huber = HuberRegressor::from(log_x2.clone(), log_y2.clone());
    let params = solve(huber);
    let alpha = 2.0 + params[0];
    (alpha, params[1])
}

struct InnerGrowth {}
