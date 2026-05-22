use std::{
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader},
    str::SplitWhitespace,
};

use anyhow::Result;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    analyses::info::FileInfo, coverage_matrix::CoverageMatrix, file_formats::FileFormatParser,
    hist::Hist,
};

pub struct VcfParser {
    filename: String,
    count_type: VcfCountType,
    /// The same as groupby_haplotypes, however,
    /// this requires the VCF to have constant ploidy for this sample
    split_haplotypes: bool,
    current_variant_id: usize,
}

impl FileFormatParser for VcfParser {
    fn generate_hist(self: Box<Self>) -> crate::hist::Hist {
        let file = File::open(&self.filename).expect("Can read file");
        let buf_reader = BufReader::new(file);
        let mut lines = buf_reader.lines().map(|l| l.expect("Failed to read line"));
        let mut header = VcfHeader::parse(&mut lines).expect("Failed to parse header");
        let paths = std::mem::take(&mut header.paths);
        let run_id = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let run_name = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let mut hist =
            Hist::from_maximum_coverage(paths.len(), self.count_type.to_string(), run_id, run_name);
        for line in lines {
            let variants = self
                .parse_variant_line_to_count(&line)
                .expect("Failed to parse variant");
            for variant in variants {
                let coverage = variant.value;
                if coverage == 0 {
                    log::info!("Variant has coverage 0: {}", variant.id);
                }
                match self.count_type {
                    VcfCountType::Variants => {
                        hist.insert_feature_of_coverage_and_length(coverage, 1);
                    }
                    VcfCountType::AltLength => {
                        hist.insert_feature_of_coverage_and_length(coverage, variant.alt_length);
                    }
                }
            }
        }
        hist
    }

    fn generate_matrix(mut self: Box<Self>) -> crate::coverage_matrix::CoverageMatrix {
        let file = File::open(&self.filename).expect("Can read file");
        let buf_reader = BufReader::new(file);
        let mut lines = buf_reader.lines().map(|l| l.expect("Failed to read line"));
        let mut header = VcfHeader::parse(&mut lines).expect("Failed to parse header");
        let paths = std::mem::take(&mut header.paths);
        let run_id = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let run_name = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let file_info = FileInfo::new("vcf");
        let mut matrix =
            CoverageMatrix::new(self.count_type.to_string(), run_id, run_name, file_info);
        let path_len = paths.len();
        let mut num_variants = 0;
        let mut num_alt_alleles = 0;
        let mut singleton_variants = 0;
        let mut first_variant = true;
        let line = lines.next().expect("Has at least one variant");
        let variants = self
            .parse_variant_line_to_allele_list(&line)
            .expect("Failed to parse variant");
        for variant in variants {
            let feature_length = match self.count_type {
                VcfCountType::Variants => 1,
                VcfCountType::AltLength => variant.alt_length,
            };
            let name = match variant.id {
                "." => self.get_variant_id(),
                id => format!("{}-{id}", self.get_variant_id()),
            };
            let pos = variant.pos;
            let feature = variant.value;
            if first_variant {
                let paths = if self.split_haplotypes {
                    paths
                        .iter()
                        .enumerate()
                        .flat_map(|(idx, name)| {
                            if feature[idx].len() == 1 {
                                vec![format!("{name}#0")]
                            } else {
                                (1..feature[idx].len() + 1)
                                    .map(|i| format!("{name}#{i}"))
                                    .collect_vec()
                            }
                        })
                        .collect()
                } else {
                    paths.clone()
                };
                matrix.set_path_names(paths);
                first_variant = false;
            }
            if feature.iter().flatten().sum::<u32>() == 1 {
                singleton_variants += 1;
            }
            matrix.insert_feature(
                name,
                feature_length,
                pos,
                feature.into_iter().flatten().collect(),
            );
            num_alt_alleles += 1;
        }
        for line in lines {
            let variants = self
                .parse_variant_line_to_allele_list(&line)
                .expect("Failed to parse variant");
            for variant in variants {
                let feature_length = match self.count_type {
                    VcfCountType::Variants => 1,
                    VcfCountType::AltLength => variant.alt_length,
                };
                let name = match variant.id {
                    "." => self.get_variant_id(),
                    id => format!("{}-{id}", self.get_variant_id()),
                };
                let pos = variant.pos;
                let feature = variant.value;
                if feature.iter().flatten().sum::<u32>() == 1 {
                    singleton_variants += 1;
                }
                matrix.insert_feature(
                    name,
                    feature_length,
                    pos,
                    feature.into_iter().flatten().collect(),
                );
                num_alt_alleles += 1;
            }
            num_variants += 1;
        }
        let mut file_info = FileInfo::new("vcf");
        file_info.add_info("Number of samples", path_len.to_string().as_str());
        file_info.add_info("Number of variants", num_variants.to_string().as_str());
        file_info.add_info(
            "Number of alt alleles",
            num_alt_alleles.to_string().as_str(),
        );
        file_info.add_info(
            "Number of singleton alleles",
            singleton_variants.to_string().as_str(),
        );
        matrix.set_file_info(file_info);
        matrix
    }
}

impl VcfParser {
    pub fn new(filename: &str, count_type: VcfCountType, split_haplotypes: bool) -> Result<Self> {
        Ok(Self {
            filename: filename.to_owned(),
            count_type,
            split_haplotypes,
            current_variant_id: 1,
        })
    }

    fn get_split_haplotype_text(&self) -> &str {
        match self.split_haplotypes {
            true => "haplotypes",
            false => "samples",
        }
    }

    fn get_variant_id(&mut self) -> String {
        let id = self.current_variant_id;
        self.current_variant_id += 1;
        id.to_string()
    }

    /// Parses a variant and returns a list noting the presence/absence of the variant for each path
    fn parse_variant_line_to_allele_list<'a>(
        &self,
        line: &'a str,
    ) -> Result<Vec<Variant<'a, Vec<Vec<u32>>>>> {
        self.parse_variant_line_base(line, |samples, num_alts, genotype_position| {
            let samples: Vec<&str> = samples.collect();
            let mut allele_lists = if self.split_haplotypes {
                vec![vec![Vec::new(); samples.len()]; num_alts]
            } else {
                vec![vec![vec![0; 1]; samples.len()]; num_alts]
            };
            samples.into_iter().enumerate().for_each(|(idx, sample)| {
                let sample_genotype = sample
                    .split(":")
                    .nth(genotype_position)
                    .expect("Sample needs a genotype");
                if !self.split_haplotypes {
                    sample_genotype
                        .split(&['|', '/'])
                        .filter_map(|x| {
                            if x == "." || x == "0" {
                                None
                            } else {
                                Some(x.parse::<usize>().expect("Allele is integer"))
                            }
                        })
                        .unique()
                        .for_each(|x| {
                            allele_lists[x - 1][idx][0] = 1;
                        });
                } else {
                    sample_genotype.split(&['|', '/']).for_each(|x| {
                        for i in 0..num_alts {
                            if x == "." {
                                allele_lists[i][idx].push(0);
                            } else {
                                let value = x.parse::<usize>().expect("Allele is integer");
                                if (i + 1) == value {
                                    allele_lists[i][idx].push(1);
                                } else {
                                    allele_lists[i][idx].push(0);
                                }
                            }
                        }
                    });
                }
            });
            allele_lists
        })
    }

    /// Parses a variant line into simple counts, how often each allele appears
    fn parse_variant_line_to_count<'a>(&self, line: &'a str) -> Result<Vec<Variant<'a, usize>>> {
        self.parse_variant_line_base(line, |samples, num_alts, genotype_position| {
            let mut counts = vec![0; num_alts];
            samples.for_each(|sample| {
                let sample_genotype = sample
                    .split(":")
                    .nth(genotype_position)
                    .expect("Sample needs a genotype");
                if !self.split_haplotypes {
                    sample_genotype
                        .split(&['|', '/'])
                        .filter_map(|x| {
                            if x == "." || x == "0" {
                                None
                            } else {
                                Some(x.parse::<usize>().expect("Allele is integer"))
                            }
                        })
                        .unique()
                        .for_each(|x| counts[x - 1] += 1);
                } else {
                    sample_genotype.split(&['|', '/']).for_each(|x| {
                        if x != "." {
                            let value = x.parse::<usize>().expect("Allele is integer");
                            if value > 0 {
                                counts[value - 1] += 1;
                            }
                        }
                    });
                }
            });
            counts
        })
    }

    /// Takes a line that contains a variant (so it does not start with #)
    /// and a closure that takes an iterator over the samples, the number of alt alleles
    /// and the position of the genotype in each sample and returns a vec with a value
    /// for each alt allele
    fn parse_variant_line_base<'a, T>(
        &self,
        line: &'a str,
        process_samples: impl FnOnce(SplitWhitespace<'a>, usize, usize) -> Vec<T>,
    ) -> Result<Vec<Variant<'a, T>>> {
        let mut fields = line.trim().split_whitespace();
        let chrom = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs chrom field!"))?;
        let pos: usize = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs pos field!"))?
            .parse()?;
        let id = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs id field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs ref field!"))?;
        let alt_lengths: Vec<usize> = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs alt field!"))?
            .split(",")
            .map(|x| if x == "." { 0 } else { x.len() })
            .collect();
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs qual field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs filter field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs info field!"))?;
        let format = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs format field!"))?;
        let genotype_position = format.split(":").position(|x| x == "GT").unwrap_or(0);
        let num_alts = alt_lengths.len();
        let values = process_samples(fields, num_alts, genotype_position);
        Ok(alt_lengths
            .into_iter()
            .zip(values.into_iter())
            .map(|(alt_length, value)| Variant {
                chrom,
                pos,
                id,
                alt_length,
                value,
            })
            .collect())
    }
}

struct Variant<'a, T> {
    #[allow(dead_code)]
    chrom: &'a str,
    pos: usize,
    id: &'a str,
    alt_length: usize,
    value: T,
}

struct VcfHeader {
    paths: Vec<String>,
}

impl VcfHeader {
    fn parse<I, S>(lines: &mut I) -> Result<Self>
    where
        I: Iterator<Item = S>,
        S: AsRef<str>,
    {
        let mut header = Self { paths: Vec::new() };

        for line in lines {
            let line = line.as_ref().trim_end();
            if line.starts_with("##") {
                continue;
            } else if line.starts_with("#CHROM") {
                let samples: Vec<&str> = line.split("\t").skip(9).collect();
                header.paths = samples.into_iter().map(|x| x.to_string()).collect();
                break;
            }
        }
        Ok(header)
    }
}

/// What is set for the length of a feature
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum VcfCountType {
    /// 1 for every feature (just count the number of variants)
    Variants,
    /// The length of the ALT allele is used
    AltLength,
}

impl Default for VcfCountType {
    fn default() -> Self {
        Self::AltLength
    }
}

impl Display for VcfCountType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Variants => write!(f, "Variant"),
            Self::AltLength => write!(f, "AltAlleleLength"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF_HEADER: &str = "##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003";

    #[test]
    fn test_vcf_header() {
        let mut header_lines = VCF_HEADER.lines();
        let vcf_header = VcfHeader::parse(&mut header_lines).unwrap();
        assert_eq!(vcf_header.paths.len(), 3);
    }

    #[test]
    fn test_vcf_parse_line_to_count() {
        let line = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
        let parser = VcfParser::new("", VcfCountType::Variants, false).unwrap();
        let mut variants = parser.parse_variant_line_to_count(line).unwrap();
        assert_eq!(variants.len(), 1);
        let variant = variants.remove(0);
        assert_eq!(variant.chrom, "20");
        assert_eq!(variant.pos, 14370);
        assert_eq!(variant.id, "rs6054257");
        assert_eq!(variant.alt_length, 1);
        assert_eq!(variant.value, 3);
    }
}
