use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::util::{CountType, Threshold};

use super::{AnalysisSection, InputRequirement};

pub struct SectionGrowth {
    count: CountType,
    sections_file: String,
    // Map from section name to group index (depends on grouping!)
    sections: HashMap<String, Vec<usize>>,
    section_order: Vec<String>,
    coverage: String,
    quorum: String,
    inner: Option<InnerGrowth>,
}

impl MatrixBasedAnalysis for SectionGrowth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        log::info!("reporting hist table");

        self.set_inner(matrix)?;
        let inner = self.inner.as_ref().unwrap();
        let mut res = String::new();

        res.push_str(&format!(
            "# {}\n",
            std::env::args().collect::<Vec<String>>().join(" ")
        ));

        res.push_str(&format!("panacus\tsection-growth\n"));
        res.push_str(&format!("count\t{}\n", self.count));

        res.push_str("coverage");
        for c in inner.thresholds.coverage.iter() {
            res.push_str(&format!("\t{}", c.get_string()));
        }
        res.push_str("\nquorum");
        for q in inner.thresholds.quorum.iter() {
            res.push_str(&format!("\t{}", q.get_string()));
        }
        res.push_str("\n");

        res.push_str(&format!("0\t{}", self.section_order[0]));
        for _ in inner.thresholds.coverage.iter() {
            res.push_str(&format!("\t{}", f64::NAN));
        }
        res.push_str("\n");
        let mut section_index = 0;
        // res.push_str(&format!("0\t{}\t0\n", inner.sections[section_index].0));
        for idx in 0..inner.growths[0].len() {
            if section_index < inner.sections.len() - 1
                && idx >= inner.sections[section_index + 1].1
            {
                section_index += 1;
            }
            res.push_str(&format!("{}\t{}", idx + 1, inner.sections[section_index].0,));
            for j in 0..inner.thresholds.coverage.len() {
                res.push_str(&format!("\t{}", inner.growths[j][idx]));
            }
            res.push_str("\n");
        }
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(matrix)?;
        Ok(vec![])
    }
}

impl SectionGrowth {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
        }
    }

    fn get_growths_for_c_q0(&self, c: Threshold, _matrix: &CoverageMatrix) -> Vec<Vec<f64>> {
        let _thresholds = ThresholdContainer {
            coverage: vec![c],
            quorum: vec![Threshold::Relative(0.0)],
        };

        unimplemented!()
    }

    fn get_growths_for_c_q(
        &self,
        _c: Threshold,
        _thresholds: ThresholdContainer,
        _matrix: &CoverageMatrix,
    ) -> Vec<Vec<f64>> {
        unimplemented!()
    }

    fn get_growths_for_coverage(
        &self,
        c: Threshold,
        q: Vec<Threshold>,
        matrix: &CoverageMatrix,
    ) -> Vec<Vec<f64>> {
        let thresholds = ThresholdContainer {
            coverage: vec![Threshold::Absolute(1); q.len()],
            quorum: q,
        };
        let (q0_thresholds, other_thresholds) = Self::split_thresholds(&thresholds);
        eprintln!(
            "q0: {}, other: {}",
            q0_thresholds.coverage.len(),
            other_thresholds.coverage.len()
        );
        // let q0_growths = self.get_growths_for_c_q0(c, gb);
        let other_growths = self.get_growths_for_c_q(c, thresholds, matrix);
        let full_growths = other_growths;
        // full_growths.extend(other_growths);

        full_growths
    }

    fn set_sections(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<Vec<(String, usize)>> {
        let sections_file = File::open(&self.sections_file)?;
        let reader = BufReader::new(sections_file);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<_> = line.trim().split_whitespace().collect();
            assert_eq!(fields.len(), 2, "Sections file can only have two columns");
            let section = fields[1].to_owned();
            if !self.section_order.contains(&section) {
                self.section_order.push(section.clone());
                self.sections.insert(section.clone(), Vec::new());
            }
            let groups = matrix.get_path_names();
            self.sections.get_mut(&section).unwrap().push(
                groups
                    .iter()
                    .position(|g| g == fields[0])
                    .expect(&format!("Group {} is not a group", fields[0])),
            );
        }
        eprintln!("section_order: {:?}", self.section_order);
        eprintln!("sections: {:?}", self.sections);

        let mut section_index = vec![];
        let mut rolling_index = 0;
        for i in 0..self.section_order.len() {
            section_index.push((self.section_order[i].clone(), rolling_index));
            rolling_index += self.sections[&self.section_order[i]].len();
        }
        Ok(section_index)
    }

    fn set_inner(&mut self, _matrix: &CoverageMatrix) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        unimplemented!();
    }

    fn split_thresholds(
        thresholds: &ThresholdContainer,
    ) -> (ThresholdContainer, ThresholdContainer) {
        let mut first_c = Vec::new();
        let mut first_q = Vec::new();
        let mut second_c = Vec::new();
        let mut second_q = Vec::new();
        for (c, q) in thresholds.coverage.iter().zip(thresholds.quorum.iter()) {
            if q.to_relative(1) == 0.0 {
                first_c.push(*c);
                first_q.push(*q);
            } else {
                second_c.push(*c);
                second_q.push(*q);
            }
        }
        (
            ThresholdContainer {
                coverage: first_c,
                quorum: first_q,
            },
            ThresholdContainer {
                coverage: second_c,
                quorum: second_q,
            },
        )
    }
}

struct InnerGrowth {
    growths: Vec<Vec<f64>>,
    sections: Vec<(String, usize)>,
    thresholds: ThresholdContainer,
}
