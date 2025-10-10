use std::cmp;
use std::collections::HashSet;
use std::io::Write;
use std::process::Command;

use once_cell::sync::Lazy;
use tempfile::NamedTempFile;

use crate::analysis_parameter::AnalysisParameter;
use crate::graph_broker::GraphBroker;
use crate::util::CountType;
use crate::{analyses::InputRequirement, html_report::ReportItem};

use super::{Analysis, AnalysisSection, ConstructibleAnalysis};

static VIRIDIS: Lazy<Vec<(f64, (f64, f64, f64))>> = Lazy::new(|| {
    vec![
        (0.0, (253.0, 231.0, 37.0)),
        (0.25, (94.0, 201.0, 98.0)),
        (0.5, (33.0, 145.0, 140.0)),
        (0.75, (59.0, 82.0, 139.0)),
        (1.0, (68.0, 1.0, 84.0)),
    ]
});

pub struct CoverageColors {
    inner: Option<InnerCoverageColors>,
}

impl ConstructibleAnalysis for CoverageColors {
    fn from_parameter(_parameter: AnalysisParameter) -> Self {
        Self { inner: None }
    }
}

impl Analysis for CoverageColors {
    fn get_type(&self) -> String {
        "CoverageColors".to_string()
    }
    fn generate_table(
        &mut self,
        gb: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        self.set_inner(gb)?;
        let mut result = String::new();
        result.push_str("Name,Color\n");
        for (id, cov) in self
            .inner
            .as_ref()
            .expect("Setting inner worked for coverage colors")
            .coverages
            .iter()
        {
            let color = convert_to_color(
                *cov,
                self.inner
                    .as_ref()
                    .expect("Setting inner worked for coverage colors")
                    .max_cov,
            );
            result.push_str(&format!("{},{}\n", id, color));
        }
        Ok(result)
    }

    fn generate_report_section(
        &mut self,
        dm: Option<&GraphBroker>,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        if dm.is_none() {
            panic!("Coverage colors needs a graph");
        }
        self.set_inner(dm)?;
        let dm = dm.unwrap();
        let table = self.generate_table(Some(dm))?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "coverage-colors-{}",
            self.get_run_id(dm)
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let output_svg = self.get_image(dm);
        let growth_tabs = vec![AnalysisSection {
            analysis: "Coverage Colors".to_string(),
            run_name: self.get_run_name(dm),
            run_id: self.get_run_id(dm),
            countable: CountType::Node.to_string(),
            table: Some(table),
            id: format!("{id_prefix}"),
            plot_downloads: vec![("svg".to_string(), "Download as svg".to_string())],
            items: match output_svg {
                Ok(svg_path) => vec![ReportItem::Svg {
                    id: format!("svg-{id_prefix}"),
                    file: svg_path,
                }],
                Err(e) => {
                    log::error!("Error running bandage: {e}");
                    Vec::new()
                }
            },
        }];
        Ok(growth_tabs)
    }

    fn get_graph_requirements(&self) -> HashSet<InputRequirement> {
        let mut req = HashSet::from([InputRequirement::AbacusByGroup(CountType::Node)]);
        req.insert(InputRequirement::Node);
        req
    }
}

impl CoverageColors {
    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }

    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-coveragecolors", gb.get_run_id())
    }

    fn set_inner(&mut self, gb: Option<&GraphBroker>) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        if gb.is_none() {
            panic!(
                "CoverageColors
             needs a graph in order to work"
            );
        }
        let gb = gb.unwrap();
        let coverages = gb.get_node_coverages();
        let node_names = gb.get_node_names();

        let coverages: Vec<(String, u64)> = coverages
            .into_iter()
            .map(|(id, cov)| (str::from_utf8(&node_names[&id]).unwrap().to_owned(), cov))
            .collect();
        let max_cov = gb.get_groups().len() as u64;

        let coverages = InnerCoverageColors { coverages, max_cov };
        self.inner = Some(coverages);
        Ok(())
    }

    fn get_image(&mut self, gb: &GraphBroker) -> anyhow::Result<String> {
        if self.inner.is_none() {
            panic!("get_image can only be called if coverage colors was already calculated");
        }

        let mut csv = NamedTempFile::with_suffix(".csv")?;
        let table = self.generate_table(Some(gb))?;
        write!(csv, "{}", table)?;
        let csv_path = csv.path().to_str().unwrap();
        let svg = NamedTempFile::with_suffix(".svg")?;
        let svg_path = svg.path().to_str().unwrap();
        let graph_path = gb.get_fname();
        Command::new("Bandage")
            .args(["image", &graph_path, svg_path, "--colors", csv_path])
            .output()?;
        let copied_path = svg_path.to_string();
        svg.keep()?;
        Ok(copied_path)
    }
}

fn convert_to_color(cov: u64, max_cov: u64) -> String {
    // Handle case of coverage == 1
    if cov == 1 {
        return "#FF0000".to_string();
    }

    if cov == 0 {
        return "#AAAAAA".to_string();
    }

    // Trim value
    let value = cmp::max(cmp::min(cov, max_cov), 2);
    let relative_value = (value - 2) as f64 / (max_cov - 2) as f64;

    for i in 0..VIRIDIS.len() - 1 {
        if relative_value >= VIRIDIS[i].0 && relative_value <= VIRIDIS[i + 1].0 {
            let factor = (relative_value - VIRIDIS[i].0) / (VIRIDIS[i + 1].0 - VIRIDIS[i].0);
            let r = (VIRIDIS[i].1 .0 + factor * (VIRIDIS[i + 1].1 .0 - VIRIDIS[i].1 .0)) as u64;
            let g = (VIRIDIS[i].1 .1 + factor * (VIRIDIS[i + 1].1 .1 - VIRIDIS[i].1 .1)) as u64;
            let b = (VIRIDIS[i].1 .2 + factor * (VIRIDIS[i + 1].1 .2 - VIRIDIS[i].1 .2)) as u64;
            return format!("#{:02X}{:02X}{:02X}", r, g, b);
        }
    }
    "#000000".to_string()
}

struct InnerCoverageColors {
    coverages: Vec<(String, u64)>,
    max_cov: u64,
}
