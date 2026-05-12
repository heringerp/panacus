use itertools::multizip;
use itertools::Itertools;

use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::{
    html_report::{AnalysisSection, Bin, ReportItem},
    util::get_default_plot_downloads,
    util::CountType,
};

pub struct NodeDistribution {
    bins: Vec<Bin>,
    min: (u32, f64),
    max: (u32, f64),
    threshold: usize,
}

impl MatrixBasedAnalysis for NodeDistribution {
    fn get_type(&self) -> String {
        "NodeDistribution".to_string()
    }

    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        if self.bins.is_empty() {
            self.set_table(matrix);
        }
        let mut result = "Bin\tCoverage\tLog-Length\tLog-Size\n".to_string();
        for (i, bin) in self.bins.iter().enumerate() {
            result.push_str(&format!("{}\t{}\t{}\t{}\n", i, bin.x, bin.y, bin.size));
        }
        Ok(result)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let table = self.generate_table(matrix)?;
        //let table = "".to_string();
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "node-dist-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let tab = vec![AnalysisSection {
            id: format!("{}-{}", id_prefix, CountType::Node.to_string()),
            analysis: "Node distribution".to_string(),
            table: Some(table),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: CountType::Node.to_string(),
            items: vec![ReportItem::Hexbin {
                id: format!("{id_prefix}-{}", CountType::Node),
                bins: self.bins.clone(),
                threshold: self.threshold,
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(tab)
    }
}

impl NodeDistribution {
    fn set_table(&mut self, matrix: &CoverageMatrix) {
        let countables = &matrix.get_feature_counts();
        let (cov_min, cov_max) = match countables.iter().minmax() {
            itertools::MinMaxResult::MinMax(min, max) => (min, max),
            _ => panic!("Node distribution needs to have at least two countables"),
        };
        let node_lens = &matrix
            .get_feature_lengths()
            .iter()
            .map(|x| (*x as f64).log10())
            .collect::<Vec<f64>>();
        let (lens_min, lens_max) = match node_lens.iter().minmax() {
            itertools::MinMaxResult::MinMax(min, max) => (min, max),
            _ => panic!("Node distribution needs to have at least two countables"),
        };
        let node_names = matrix.get_feature_names();
        let points: Vec<(String, u32, f64)> = multizip((
            (0..node_lens.len())
                .into_iter()
                .map(|id| (&node_names[id]).to_string())
                .collect::<Vec<String>>(),
            countables.into_iter().copied().map(|x| x as u32),
            node_lens.into_iter().copied(),
        ))
        .collect();
        let bins = Bin::hexbin(&points, 15, 9);
        self.bins = bins;
        self.min = (*cov_min as u32, *lens_min);
        self.max = (*cov_max as u32, *lens_max);
    }
}
