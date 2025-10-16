use core::panic;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::analysis_parameter::AnalysisParameter;
use crate::graph_broker::GraphBroker;
use crate::util::CountType;

use super::{Analysis, AnalysisSection, ConstructibleAnalysis, InputRequirement};

pub struct SectionGrowth {
    count: CountType,
    sections_file: String,
    // Map from section name to group index (depends on grouping!)
    sections: HashMap<String, Vec<usize>>,
    section_order: Vec<String>,
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
        res.push_str(&format!("coverage\t1\n"));
        res.push_str(&format!("quorum\t0\n"));

        let mut section_index = 0;
        res.push_str(&format!("0\t{}\t0\n", inner.sections[section_index].0));
        for (idx, growth) in inner.growth.iter().enumerate() {
            if section_index < inner.sections.len() - 1
                && idx >= inner.sections[section_index + 1].1
            {
                section_index += 1;
            }
            res.push_str(&format!(
                "{}\t{}\t{}\n",
                idx + 1,
                inner.sections[section_index].0,
                growth
            ));
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
        SectionGrowth {
            count: match &parameter {
                AnalysisParameter::SectionGrowth { count_type, .. } => *count_type,
                _ => panic!("Similarity analysis needs similarity parameter"),
            },
            sections_file: match &parameter {
                AnalysisParameter::SectionGrowth { sections, .. } => sections.to_owned(),
                _ => panic!("Similarity analysis needs similarity parameter"),
            },
            sections: HashMap::new(),
            section_order: Vec::new(),
            inner: None,
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

    fn set_inner(&mut self, gb: Option<&GraphBroker>) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        if gb.is_none() {
            panic!("Section Growth needs a graph");
        }
        let gb = gb.unwrap();
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

        let first_group = &self.sections[&self.section_order[0]];
        let (abacus, non_zero_nodes) = gb.get_abacus_by_group().to_abacus_all_nodes(first_group);
        let hist = gb.get_hist_from_abacus(&abacus);
        eprintln!("hist: {:?}", hist);
        let mut full_growth = hist.calc_growth(
            &crate::util::Threshold::Absolute(1),
            &crate::util::Threshold::Relative(0.0),
        );
        eprintln!("growth: {:?}", full_growth);
        let mut section_index = vec![(self.section_order[0].clone(), 0)];
        let mut current_index = first_group.len();
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
            let growth = hist.calc_growth(
                &crate::util::Threshold::Absolute(1),
                &crate::util::Threshold::Relative(0.0),
            );
            let growth: Vec<f64> = growth
                .into_iter()
                .map(|x| x + full_growth[full_growth.len() - 1])
                .collect();
            full_growth.extend(growth.into_iter());
            already_seen_nodes.extend(non_zero_nodes.into_iter());
            section_index.push((section.clone(), current_index));
            current_index += paths.len();
        }
        self.inner = Some(InnerGrowth {
            growth: full_growth,
            sections: section_index,
        });

        Ok(())
    }
}

struct InnerGrowth {
    growth: Vec<f64>,
    sections: Vec<(String, usize)>,
}
