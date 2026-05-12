use anyhow::Error;
use itertools::Itertools;

use crate::{
    coverage_matrix::CoverageMatrix,
    file_formats::{
        gfa_parser::{
            graph::GraphStorage,
            util::{parse_path_seq_update_tables, parse_walk_seq_update_tables, update_tables},
        },
        FileFormatParser,
    },
    hist::Hist,
    util::{CountType, GroupSize, ItemIdSize},
};

use std::io::BufRead;
use std::time::Instant;
use std::{collections::HashMap, str};

use abacus::AbacusByTotal;
use util::{
    parse_path_identifier, parse_path_seq_to_item_vec, parse_walk_identifier,
    parse_walk_seq_to_item_vec, update_tables_edgecount,
};

use crate::io::bufreader_from_compressed_gfa;
use crate::util::{intersects, is_contained, ActiveTable, IntervalContainer, ItemTable};

mod abacus;
mod graph;
mod hist;
mod sparse_matrix;
mod util;

pub use abacus::GraphMask;
pub use abacus::GraphMaskParameters;
pub use graph::Edge;
pub use graph::ItemId;
pub use graph::Orientation;
pub use graph::PathSegment;
pub use hist::Hist3D;
pub use hist::ThresholdContainer;
pub use sparse_matrix::SparseMatrix;

pub struct GfaParser {
    // Inputs
    filename: String,
    count_type: CountType,
    graph_mask_parameters: GraphMaskParameters,
    is_nice: bool,

    // Generated
    graph_storage: GraphStorage,
    graph_mask: GraphMask,
}

impl FileFormatParser for GfaParser {
    fn generate_hist(self: Box<Self>) -> Hist {
        let number_of_groups = self.graph_mask.count_groups();
        let mut hist = Hist::from_maximum_coverage(
            number_of_groups,
            self.count_type.to_string(),
            self.get_run_id(),
            self.get_run_name(),
        );
        let (abacus, _, _) = self.get_abacus_by_total();
        for feature_coverage in abacus.countable {
            hist.insert_feature_of_coverage(feature_coverage as usize);
        }
        hist
    }

    fn generate_matrix(self: Box<Self>) -> CoverageMatrix {
        let (item_table, path_names, feature_lengths) = self
            .get_cleaned_item_table(
                &self.graph_mask,
                &self.graph_storage,
                self.count_type,
                &Vec::new(),
            )
            .expect("Can parse GFA file");
        let mut matrix = CoverageMatrix::new(
            self.count_type.to_string(),
            self.get_run_id(),
            self.get_run_name(),
        );
        let number_of_features = feature_lengths.len();
        matrix.insert_item_table(
            path_names,
            feature_lengths,
            vec![0; number_of_features],
            item_table,
        );
        matrix
    }

    fn get_file_info(&self) {
        unimplemented!()
    }
}

impl GfaParser {
    pub fn new(
        filename: &str,
        count_type: CountType,
        graph_mask_parameters: GraphMaskParameters,
        is_nice: bool,
    ) -> Result<Self, Error> {
        let graph_storage = GraphStorage::from_gfa(filename, is_nice, count_type);
        let graph_mask = GraphMask::from_datamgr(&graph_mask_parameters, &graph_storage)?;
        Ok(Self {
            filename: filename.to_owned(),
            count_type,
            graph_mask_parameters,
            is_nice,

            graph_storage,
            graph_mask,
        })
    }

    fn get_run_name(&self) -> String {
        format!(
            "{} {} {}",
            self.filename,
            self.graph_mask_parameters.positive_list,
            self.graph_mask_parameters.groupby
        )
    }

    fn get_run_id(&self) -> String {
        format!(
            "{}-{}-{}",
            self.filename,
            self.graph_mask_parameters.positive_list,
            self.graph_mask_parameters.groupby
        )
    }

    pub fn get_abacus_by_total(
        &self,
    ) -> (
        AbacusByTotal,
        HashMap<PathSegment, (u32, u32)>,
        HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    ) {
        let (item_table, exclude_table, subset_covered_bps, path_lens, collected_paths) =
            self.parse_paths_walks(&self.count_type, &Vec::new());
        let abacus = AbacusByTotal::item_table_to_abacus(
            &self.graph_mask,
            &self.graph_storage,
            self.count_type,
            item_table,
            exclude_table,
            subset_covered_bps,
        );
        (abacus, path_lens, collected_paths)
    }

    fn get_cleaned_item_table(
        &self,
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        count: CountType,
        paths_to_collect: &Vec<PathSegment>,
    ) -> anyhow::Result<(ItemTable, Vec<String>, Vec<usize>)> {
        log::info!("parsing path + walk sequences");
        let (item_table, exclude_table, subset_covered_bps, _paths_len, _collected_paths) =
            self.parse_paths_walks(&count, paths_to_collect);

        let mut path_order: Vec<(ItemIdSize, GroupSize)> = Vec::new();
        let mut groups: Vec<String> = Vec::new();

        for (path_id, group_id) in graph_mask.get_path_order(&graph_storage.path_segments) {
            log::debug!(
                "processing path {} (group {})",
                &graph_storage.path_segments[path_id as usize],
                group_id
            );
            if groups.is_empty() || groups.last().unwrap() != group_id {
                groups.push(group_id.to_string());
            }
            //if groups.len() > 65534 {
            //    panic!("data has more than 65534 path groups, but command is not supported for more than 65534");
            //}
            path_order.push((path_id, (groups.len() - 1) as GroupSize));
        }

        println!("EXCLUDE_TABLE: {:?}", exclude_table);
        println!("SUBSET_COVERED_BPS: {:?}", subset_covered_bps);

        let path_names = graph_mask
            .groups
            .iter()
            .map(|x| x.0.to_string())
            .collect_vec();
        let feature_lengths = match count {
            CountType::Node => vec![1; graph_storage.node_count],
            CountType::Edge => vec![1; graph_storage.edge_count],
            CountType::Bp => graph_storage
                .node_lens
                .iter()
                .map(|x| *x as usize)
                .collect_vec(),
        };

        Ok((item_table, path_names, feature_lengths))
    }

    fn parse_paths_walks(
        &self,
        count_type: &CountType,
        paths_to_collect: &Vec<PathSegment>,
    ) -> (
        ItemTable,
        Option<ActiveTable>,
        Option<IntervalContainer>,
        HashMap<PathSegment, (u32, u32)>,
        HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    ) {
        log::info!("parsing path + walk sequences");
        log::info!("collecting: {:?}", paths_to_collect);
        let mut data = bufreader_from_compressed_gfa(&self.filename);
        let graph_storage = &self.graph_storage;
        let graph_mask = &self.graph_mask;
        let mut item_table = ItemTable::new(graph_storage.path_segments.len());

        let (mut subset_covered_bps, mut exclude_table, include_map, exclude_map) =
            graph_mask.load_optional_subsetting(graph_storage, count_type);

        let mut collected_paths: HashMap<PathSegment, Vec<(ItemId, Orientation)>> = HashMap::new();

        let mut num_path = 0;
        let complete: Vec<(usize, usize)> = vec![(0, usize::MAX)];
        let mut paths_len: HashMap<PathSegment, (u32, u32)> = HashMap::new();

        let mut buf = vec![];
        let timer = Instant::now();
        while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
            if buf[0] == b'P' || buf[0] == b'W' {
                let (path_seg, buf_path_seg) = match buf[0] {
                    b'P' => parse_path_identifier(&buf),
                    b'W' => parse_walk_identifier(&buf),
                    _ => unreachable!(),
                };

                log::debug!("processing path {}", &path_seg);

                let include_coords =
                    get_include_coords(graph_mask, &complete, &include_map, &path_seg);
                let exclude_coords = get_exclude_maps(graph_mask, &exclude_map, &path_seg);

                let (start, end) = path_seg.coords().unwrap_or((0, usize::MAX));

                // do not process the path sequence if path is neither part of subset nor exclude
                if graph_mask.include_coords.is_some()
                    && !intersects(include_coords, &(start, end))
                    && !intersects(exclude_coords, &(start, end))
                {
                    log::debug!("path {} does not intersect with subset coordinates {:?} nor with exclude coordinates {:?} and therefore is skipped from processing", &path_seg, &include_coords, &exclude_coords);
                    skip_path(&mut item_table, &mut num_path, &mut buf);
                    continue;
                }

                if *count_type != CountType::Edge
                    && (graph_mask.include_coords.is_none()
                        || is_contained(include_coords, &(start, end)))
                    && (graph_mask.exclude_coords.is_none()
                        || is_contained(exclude_coords, &(start, end)))
                {
                    log::debug!("path {} is fully contained within subset coordinates {:?} and is eligible for full parallel processing", path_seg, include_coords);
                    let mut none = None;
                    let ex: &mut Option<ActiveTable> = if exclude_coords.is_empty() {
                        &mut none
                    } else {
                        &mut exclude_table
                    };
                    let (num_added_nodes, bp_len) = match buf[0] {
                        b'P' => parse_path_seq_update_tables(
                            buf_path_seg,
                            graph_storage,
                            &mut item_table,
                            ex.as_mut(),
                            num_path,
                        ),
                        b'W' => parse_walk_seq_update_tables(
                            buf_path_seg,
                            graph_storage,
                            &mut item_table,
                            ex.as_mut(),
                            num_path,
                        ),
                        _ => unreachable!(),
                    };
                    paths_len.insert(path_seg.clone(), (num_added_nodes, bp_len));
                } else {
                    let sids = match buf[0] {
                        b'P' => parse_path_seq_to_item_vec(buf_path_seg, graph_storage),
                        b'W' => parse_walk_seq_to_item_vec(buf_path_seg, graph_storage),
                        _ => unreachable!(),
                    };
                    if paths_to_collect.iter().any(|p| path_seg.is_part_of(p)) {
                        collected_paths.insert(path_seg.clone(), sids.clone());
                    }
                    match count_type {
                        CountType::Node | CountType::Bp => {
                            let (node_len, bp_len) = update_tables(
                                &mut item_table,
                                &mut subset_covered_bps.as_mut(),
                                &mut exclude_table.as_mut(),
                                num_path,
                                graph_storage,
                                sids,
                                include_coords,
                                exclude_coords,
                                start,
                            );
                            paths_len.insert(path_seg.clone(), (node_len as u32, bp_len as u32));
                        }
                        CountType::Edge => update_tables_edgecount(
                            &mut item_table,
                            &mut exclude_table.as_mut(),
                            num_path,
                            graph_storage,
                            sids,
                            include_coords,
                            exclude_coords,
                            start,
                        ),
                    };
                }
                num_path += 1;
            }
            buf.clear();
        }
        let duration = timer.elapsed();
        log::info!(
            "func done; count: {:?}; time elapsed: {:?}",
            count_type,
            duration
        );

        (
            item_table,
            exclude_table,
            subset_covered_bps,
            paths_len,
            collected_paths,
        )
    }
}

fn get_exclude_maps<'a>(
    graph_mask: &GraphMask,
    exclude_map: &'a HashMap<String, Vec<(usize, usize)>>,
    path_seg: &PathSegment,
) -> &'a [(usize, usize)] {
    let exclude_coords = if graph_mask.exclude_coords.is_none() {
        &[]
    } else {
        match exclude_map.get(&path_seg.id()) {
            None => &[],
            Some(coords) => {
                log::debug!(
                    "found exclude coords {:?} for path segment {}",
                    &coords[..],
                    &path_seg.id()
                );
                &coords[..]
            }
        }
    };
    exclude_coords
}

fn get_include_coords<'a>(
    graph_mask: &GraphMask,
    complete: &'a Vec<(usize, usize)>,
    include_map: &'a HashMap<String, Vec<(usize, usize)>>,
    path_seg: &PathSegment,
) -> &'a [(usize, usize)] {
    if graph_mask.include_coords.is_none() {
        complete
    } else {
        match include_map.get(&path_seg.id()) {
            None => &[],
            Some(coords) => {
                log::debug!(
                    "found include coords {:?} for path segment {}",
                    &coords[..],
                    path_seg.id()
                );
                &coords[..]
            }
        }
    }
}

fn get_count_indices(count_types: &Vec<CountType>) -> HashMap<CountType, Vec<usize>> {
    let mut indices: HashMap<CountType, Vec<usize>> = HashMap::new();
    for (i, count_type) in count_types
        .iter()
        .map(|c| match c {
            CountType::Node => &CountType::Bp,
            count => count,
        })
        .enumerate()
    {
        if let Some(entry) = indices.get_mut(count_type) {
            (*entry).push(i);
        } else {
            indices.insert(*count_type, vec![i]);
        }
    }
    indices
}

fn skip_path(item_table: &mut ItemTable, num_path: &mut usize, buf: &mut Vec<u8>) {
    item_table.id_prefsum[*num_path + 1] += item_table.id_prefsum[*num_path];
    *num_path += 1;
    buf.clear();
}
