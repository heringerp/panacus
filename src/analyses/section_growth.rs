use core::panic;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::analysis_parameter::AnalysisParameter;
use crate::graph_broker::{GraphBroker, ThresholdContainer};
use crate::util::{CountType, Threshold};

use super::{Analysis, AnalysisSection, ConstructibleAnalysis, InputRequirement};

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

impl Analysis for SectionGrowth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(
        &mut self,
        dm: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        log::info!("reporting hist table");

        self.set_inner(dm)?;
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
        dm: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(dm)?;
        Ok(vec![])
    }

    fn get_graph_requirements(&self) -> HashSet<super::InputRequirement> {
        let mut req = HashSet::from([InputRequirement::AbacusByGroup(self.count)]);
        req.extend(Self::count_to_input_req(self.count));
        req
    }
}

impl ConstructibleAnalysis for SectionGrowth {
    fn from_parameter(parameter: AnalysisParameter) -> Self {
        if let AnalysisParameter::SectionGrowth {
            count_type,
            sections,
            coverage,
            quorum,
        } = &parameter
        {
            SectionGrowth {
                count: *count_type,
                sections_file: sections.to_owned(),
                sections: HashMap::new(),
                section_order: Vec::new(),
                quorum: quorum.to_owned().unwrap_or("0".to_string()),
                coverage: coverage.to_owned().unwrap_or("1".to_string()),
                inner: None,
            }
        } else {
            panic!("SectionGrowth analysis needs section growth parameter");
        }
    }
}

impl SectionGrowth {
    fn count_to_input_req(count: CountType) -> HashSet<InputRequirement> {
        match count {
            CountType::Bp => HashSet::from([InputRequirement::Bp]),
            CountType::Node => HashSet::from([InputRequirement::Node]),
            CountType::Edge => HashSet::from([InputRequirement::Edge]),
            CountType::All => HashSet::from([
                InputRequirement::Bp,
                InputRequirement::Node,
                InputRequirement::Edge,
            ]),
        }
    }

    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }
    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-growth", gb.get_run_id())
    }

    fn get_growths_for_c_q0(&self, c: Threshold, gb: &GraphBroker) -> Vec<Vec<f64>> {
        let thresholds = ThresholdContainer {
            coverage: vec![c],
            quorum: vec![Threshold::Relative(0.0)],
        };

        let first_group = &self.sections[&self.section_order[0]];
        let (abacus, non_zero_nodes) = gb.get_abacus_by_group().to_abacus_all_nodes(first_group);
        let hist = gb.get_hist_from_abacus(&abacus);
        eprintln!("hist: {:?}", hist);
        let mut full_growths = hist.calc_all_growths(&thresholds, false);
        eprintln!("growth: {:?}", full_growths);
        let mut already_seen_nodes: HashSet<usize> = non_zero_nodes.into_iter().collect();
        for section in self.section_order.iter().skip(1) {
            log::info!("Handling section {}", section);
            let paths = &self.sections[section];
            let nodes: HashSet<usize> = gb
                .get_abacus_by_group()
                .get_nodes_of_paths(paths)
                .into_iter()
                .collect();
            log::info!("Gotten {} nodes of section", nodes.len());
            let nodes: Vec<usize> = nodes
                .difference(&already_seen_nodes)
                .into_iter()
                .copied()
                .collect();
            log::info!("Working with {} nodes after difference", nodes.len());
            let (abacus, non_zero_nodes) = gb.get_abacus_by_group().to_abacus(&nodes, paths);
            log::info!("Calculating another hist");
            let hist = gb.get_hist_from_abacus(&abacus);
            log::info!("Calculating another growth");
            let growths = hist.calc_all_growths(&thresholds, false);
            let growths: Vec<Vec<f64>> = growths
                .into_iter()
                .enumerate()
                .map(|(idx, growth)| {
                    growth
                        .into_iter()
                        .map(|x| x + full_growths[idx][full_growths[idx].len() - 1])
                        .collect()
                })
                .collect();
            full_growths = full_growths
                .into_iter()
                .zip(growths.into_iter())
                .map(|(mut full_growth, growth)| {
                    full_growth.extend(growth.into_iter());
                    full_growth
                })
                .collect();
            already_seen_nodes.extend(non_zero_nodes.into_iter());
        }
        full_growths
    }

    fn get_growths_for_c_q(
        &self,
        c: Threshold,
        thresholds: ThresholdContainer,
        gb: &GraphBroker,
    ) -> Vec<Vec<f64>> {
        let first_group = &self.sections[&self.section_order[0]];
        let (abacus, _) = gb
            .get_abacus_by_group()
            .to_abacus_all_nodes_with_coverage(first_group, &c);
        let hist = gb.get_hist_from_abacus(&abacus);
        let mut full_growths = hist.calc_all_growths(&thresholds, false);
        let mut path_collection: Vec<usize> = first_group.clone();

        for section in self.section_order.iter().skip(1) {
            log::info!("Handling quorum section {}", section);
            let paths = &self.sections[section];
            let (current_abacus, _) = gb
                .get_abacus_by_group()
                .to_abacus_all_nodes_with_coverage(paths, &c);
            let (previous_abacus, _) = gb
                .get_abacus_by_group()
                .to_abacus_all_nodes_with_coverage(&path_collection, &c);
            let hist = gb.get_hist3d_from_abaci(&previous_abacus, &current_abacus);
            let growths = hist.calc_all_growths(&thresholds, false);
            full_growths = full_growths
                .into_iter()
                .zip(growths.into_iter())
                .map(|(mut full_growth, growth)| {
                    full_growth.extend(growth.into_iter());
                    full_growth
                })
                .collect();
            path_collection.extend(paths.into_iter());
        }
        full_growths
    }

    fn get_growths_for_coverage(
        &self,
        c: Threshold,
        q: Vec<Threshold>,
        gb: &GraphBroker,
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
        let other_growths = self.get_growths_for_c_q(c, thresholds, gb);
        let full_growths = other_growths;
        // full_growths.extend(other_growths);

        full_growths
    }

    fn set_sections(&mut self, gb: &GraphBroker) -> anyhow::Result<Vec<(String, usize)>> {
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
            let groups = gb.get_abacus_by_group().get_groups();
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

    fn set_inner(&mut self, gb: Option<&GraphBroker>) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        if gb.is_none() {
            panic!("Section Growth needs a graph");
        }
        let gb = gb.unwrap();
        let thresholds = ThresholdContainer::parse_params(&self.quorum, &self.coverage)?;

        let section_index = self.set_sections(gb)?;

        let coverage_quorums: HashMap<usize, Vec<_>> = thresholds
            .coverage
            .iter()
            .zip(thresholds.quorum.iter())
            .fold(HashMap::new(), |mut acc, (c, q)| match c {
                Threshold::Absolute(c) => {
                    acc.entry(*c).or_default().push(*q);
                    acc
                }
                _ => panic!("Coverage can not be relative for section growth"),
            });
        let mut full_growths = Vec::new();
        let mut thresholds_sorted = Vec::new();
        for (c, qs) in coverage_quorums.into_iter() {
            let c = Threshold::Absolute(c);
            let new_thresholds: Vec<(Threshold, Threshold)> = qs.iter().map(|q| (c, *q)).collect();
            thresholds_sorted.extend(new_thresholds.into_iter());
            full_growths.extend(self.get_growths_for_coverage(c, qs, gb).into_iter());
        }
        let (coverage_sorted, quorum_sorted): (Vec<_>, Vec<_>) =
            thresholds_sorted.into_iter().unzip();
        let thresholds_sorted = ThresholdContainer {
            coverage: coverage_sorted,
            quorum: quorum_sorted,
        };

        self.inner = Some(InnerGrowth {
            growths: full_growths,
            sections: section_index,
            thresholds: thresholds_sorted,
        });

        Ok(())
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
