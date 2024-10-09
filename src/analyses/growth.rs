use std::io::Write;
use std::{
    collections::HashSet,
    fs,
    io::{BufReader, BufWriter, Error},
};

use clap::{arg, value_parser, Arg, Command};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::data_manager::Hist;
use crate::{clap_enum_variants, io::OutputFormat};
use crate::{
    data_manager::{HistAuxilliary, ViewParams},
    io::{parse_hists, write_table},
    util::CountType,
};

use super::{Analysis, ReportSection};

pub struct Growth {
    growths: Vec<(CountType, Vec<Vec<f64>>)>,
    comments: Vec<Vec<u8>>,
    hist_aux: HistAuxilliary,
    hists: Vec<Hist>,
}

impl Analysis for Growth {
    fn build(
        _dm: &crate::data_manager::DataManager,
        matches: &clap::ArgMatches,
    ) -> Result<Box<Self>, Error> {
        let matches = matches.subcommand_matches("growth").unwrap();
        let coverage = matches.get_one::<String>("coverage").cloned().unwrap();
        let quorum = matches.get_one::<String>("quorum").cloned().unwrap();
        let hist_aux = HistAuxilliary::parse_params(&quorum, &coverage)?;
        let hist_file = matches
            .get_one::<String>("hist_file")
            .cloned()
            .unwrap_or_default();
        log::info!("loading coverage histogram from {}", hist_file);
        let mut data = BufReader::new(fs::File::open(&hist_file)?);
        let (coverages, comments) = parse_hists(&mut data)?;
        let hists: Vec<Hist> = coverages
            .into_iter()
            .map(|(count, coverage)| Hist { count, coverage })
            .collect();
        let growths: Vec<(CountType, Vec<Vec<f64>>)> = hists
            .par_iter()
            .map(|h| (h.count, h.calc_all_growths(&hist_aux)))
            .collect();
        Ok(Box::new(Self {
            growths,
            comments,
            hist_aux,
            hists,
        }))
    }

    fn write_table<W: Write>(
        &mut self,
        _dm: &crate::data_manager::DataManager,
        out: &mut BufWriter<W>,
    ) -> Result<(), Error> {
        log::info!("reporting hist table");
        for c in &self.comments {
            out.write_all(&c[..])?;
            out.write_all(b"\n")?;
        }
        writeln!(
            out,
            "# {}",
            std::env::args().collect::<Vec<String>>().join(" ")
        )?;

        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            "coverage".to_string(),
            "quorum".to_string(),
        ]];
        let mut output_columns: Vec<Vec<f64>> = Vec::new();

        for h in &self.hists {
            output_columns.push(h.coverage.iter().map(|x| *x as f64).collect());
            header_cols.push(vec![
                "hist".to_string(),
                h.count.to_string(),
                String::new(),
                String::new(),
            ])
        }

        for (count, g) in &self.growths {
            output_columns.extend(g.clone());
            let m = self.hist_aux.coverage.len();
            header_cols.extend(
                std::iter::repeat("growth")
                    .take(m)
                    .zip(std::iter::repeat(count).take(m))
                    .zip(self.hist_aux.coverage.iter())
                    .zip(&self.hist_aux.quorum)
                    .map(|(((p, t), c), q)| {
                        vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                    }),
            );
        }
        write_table(&header_cols, &output_columns, out)
    }

    fn generate_report_section(
        &mut self,
        _dm: &crate::data_manager::DataManager,
    ) -> super::ReportSection {
        ReportSection {}
    }

    fn get_subcommand() -> Command {
        Command::new("growth")
            .about("Calculate growth curve from coverage histogram")
            .args(&[
                arg!(hist_file: <HIST_FILE> "Coverage histogram as tab-separated value (tsv) file"),
                arg!(-a --hist "Also include histogram in output"),
                Arg::new("output_format").help("Choose output format: table (tab-separated-values) or html report").short('o').long("output-format")
                    .default_value("table").value_parser(clap_enum_variants!(OutputFormat)).ignore_case(true),
                Arg::new("coverage").help("Ignore all countables with a coverage lower than the specified threshold. The coverage of a countable corresponds to the number of path/walk that contain it. Repeated appearances of a countable in the same path/walk are counted as one. You can pass a comma-separated list of coverage thresholds, each one will produce a separated growth curve (e.g., --coverage 2,3). Use --quorum to set a threshold in conjunction with each coverage (e.g., --quorum 0.5,0.9)")
                    .short('l').long("coverage").default_value("1"),
                Arg::new("quorum").help("Unlike the --coverage parameter, which specifies a minimum constant number of paths for all growth point m (1 <= m <= num_paths), --quorum adjust the threshold based on m. At each m, a countable is counted in the average growth if the countable is contained in at least floor(m*quorum) paths. Example: A quorum of 0.9 requires a countable to be in 90% of paths for each subset size m. At m=10, it must appear in at least 9 paths. At m=100, it must appear in at least 90 paths. A quorum of 1 (100%) requires presence in all paths of the subset, corresponding to the core. Default: 0, a countable counts if it is present in any path at each growth point. Specify multiple quorum values with a comma-separated list (e.g., --quorum 0.5,0.9). Use --coverage to set static path thresholds in conjunction with variable quorum percentages (e.g., --coverage 5,10).")
                    .short('q').long("quorum").default_value("0"),
                Arg::new("threads").short('t').long("threads").help("").default_value("0").value_parser(value_parser!(usize)),
            ])
    }

    fn get_input_requirements(
        _matches: &clap::ArgMatches,
    ) -> Option<(HashSet<super::InputRequirement>, ViewParams, String)> {
        None
    }
}