use crate::{clap_enum_variants, util::CountType};
use clap::{arg, Arg, ArgMatches, Command};

use crate::analysis_parameter::{AnalysisParameter, AnalysisRun, Grouping};

pub fn get_subcommand() -> Command {
    Command::new("section-growth")
        .about("Calculate growth curve for sections of paths")
        .args(&[
            arg!(file: <FILE> "A graph in GFA1 format, accepts also compressed (.gz) file"),
            arg!(section_file: <FILE> "A TSV file containing paths in one column and the corresponding sections in the other"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file) (ONLY IN GFA MODE)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list (ONLY IN GFA MODE)"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file (ONLY IN GFA MODE)"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype (ONLY IN GFA MODE)"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample (ONLY IN GFA MODE)"),
            Arg::new("count").help("Graph quantity to be counted").default_value("node").ignore_case(true).short('c').long("count").value_parser(clap_enum_variants!(CountType)),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<Result<Vec<AnalysisRun>, anyhow::Error>> {
    if let Some(args) = args.subcommand_matches("section-growth") {
        let section_file = args
            .get_one::<String>("section_file")
            .expect("SectionGrowth subcommand has section file")
            .to_owned();
        let count = args
            .get_one::<CountType>("count")
            .expect("SectionGrowth subcommand has count type")
            .to_owned();
        let graph = args
            .get_one::<String>("file")
            .expect("SectionGrowth subcommand has gfa file")
            .to_owned();
        let subset = args
            .get_one::<String>("subset")
            .cloned()
            .unwrap_or_default();
        let exclude = args
            .get_one::<String>("exclude")
            .cloned()
            .unwrap_or_default();
        let grouping = args.get_one::<String>("groupby").cloned();
        let grouping = if args.get_flag("groupby-sample") {
            Some(Grouping::Sample)
        } else if args.get_flag("groupby-haplotype") {
            Some(Grouping::Haplotype)
        } else {
            grouping.map(|g| Grouping::Custom(g))
        };
        Some(Ok(vec![AnalysisRun::new(
            graph,
            None,
            subset,
            exclude,
            grouping,
            false,
            vec![AnalysisParameter::SectionGrowth {
                sections: section_file,
                count_type: count,
            }],
        )]))
    } else {
        None
    }
}
