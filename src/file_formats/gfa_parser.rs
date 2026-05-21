use anyhow::Error;
use itertools::Itertools;

use crate::{
    analyses::info::FileInfo,
    coverage_matrix::CoverageMatrix,
    file_formats::{
        gfa_parser::{
            graph::GraphStorage,
            util::{parse_path_seq_update_tables, parse_walk_seq_update_tables, update_tables},
        },
        FileFormatParser,
    },
    hist::Hist,
    util::{
        averageu32, median_already_sorted, n50_already_sorted, CountType, GroupSize, ItemIdSize,
    },
};

use std::time::Instant;
use std::{collections::HashMap, str};
use std::{collections::HashSet, io::BufRead};

pub use abacus::AbacusByTotal;
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
pub use hist::choose;
pub use hist::Hist3D;
pub use hist::ThresholdContainer;
pub use sparse_matrix::SparseMatrix;

pub struct GfaParser {
    // Inputs
    filename: String,
    count_type: CountType,
    graph_mask_parameters: GraphMaskParameters,

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
        for (idx, feature_coverage) in abacus.countable.iter().enumerate().skip(1) {
            if self.count_type == CountType::Bp {
                let length = self.graph_storage.node_lens[idx];
                hist.insert_feature_of_coverage_and_length(
                    *feature_coverage as usize,
                    length as usize,
                );
            } else {
                hist.insert_feature_of_coverage(*feature_coverage as usize);
            }
        }
        hist
    }

    fn generate_matrix(self: Box<Self>) -> CoverageMatrix {
        let (item_table, path_names, feature_lengths, feature_names, path_lengths) = self
            .get_cleaned_item_table(
                &self.graph_mask,
                &self.graph_storage,
                self.count_type,
                &Vec::new(),
            )
            .expect("Can parse GFA file");
        let file_info = self.get_file_info_value(path_lengths);
        let mut matrix = CoverageMatrix::new(
            self.count_type.to_string(),
            self.get_run_id(),
            self.get_run_name(),
            file_info,
        );
        let number_of_features = feature_lengths.len();
        matrix.insert_item_table(
            path_names,
            feature_lengths,
            vec![0; number_of_features],
            feature_names,
            item_table,
        );
        matrix
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
            graph_storage,
            graph_mask,
        })
    }

    fn get_file_info_value(&self, paths_len: HashMap<PathSegment, (u32, u32)>) -> FileInfo {
        // Created here to satisfy the borrow checker
        let empty_degree = vec![0, 0];
        let degree = self.graph_storage.degree.as_ref().unwrap_or(&empty_degree);
        let mut node_lens_sorted = self.graph_storage.node_lens[1..].to_vec();
        node_lens_sorted.sort_by(|a, b| b.cmp(a)); // decreasing, for N50
        let mut components = connected_components(
            self.graph_storage
                .edge2id
                .as_ref()
                .unwrap_or(&HashMap::new()),
            &self.graph_storage.get_nodes(),
        );
        components.sort();

        let paths_bp_len: Vec<_> = paths_len.values().map(|x| x.1).collect();
        let paths_len: Vec<_> = paths_len.values().map(|x| x.0).collect();

        // group_count: gb.get_group_count(),
        let mut file_info = FileInfo::new("gfa");
        file_info.add_info(
            "Number of nodes",
            self.graph_storage.node_count.to_string().as_str(),
        );
        file_info.add_info(
            "Number of edges",
            self.graph_storage.edge_count.to_string().as_str(),
        );
        file_info.add_info(
            "Average degree",
            averageu32(&degree[1..]).to_string().as_str(),
        );
        file_info.add_info(
            "Max degree",
            degree[1..].iter().max().unwrap().to_string().as_str(),
        );
        file_info.add_info(
            "Min degree",
            degree[1..].iter().min().unwrap().to_string().as_str(),
        );
        file_info.add_info(
            "Number of zero-degree nodes",
            degree[1..]
                .iter()
                .filter(|&x| *x == 0)
                .count()
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Number of connected components",
            (components.len() as u32).to_string().as_str(),
        );
        file_info.add_info(
            "Largest component size",
            (*components.iter().max().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Smallest component size",
            (*components.iter().min().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Median component size",
            (median_already_sorted(&components)).to_string().as_str(),
        );
        file_info.add_info(
            "Longest node",
            (*node_lens_sorted.iter().max().unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Shortest node",
            (*node_lens_sorted.iter().min().unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Average node size",
            (averageu32(&node_lens_sorted)).to_string().as_str(),
        );
        file_info.add_info(
            "Median node size",
            (median_already_sorted(&node_lens_sorted))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "N50 nodes",
            (n50_already_sorted(&node_lens_sorted).unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Number of basepairs",
            (node_lens_sorted.iter().sum::<u32>()).to_string().as_str(),
        );
        file_info.add_info(
            "Number of groups",
            (self.graph_mask.count_groups()).to_string().as_str(),
        );
        file_info.add_info(
            "Number of paths",
            (self.graph_storage.path_segments.len())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Longest path in nodes",
            (*paths_len.iter().max().unwrap()).to_string().as_str(),
        );
        file_info.add_info(
            "Shortest path in nodes",
            (*paths_len.iter().min().unwrap()).to_string().as_str(),
        );
        file_info.add_info(
            "Average path length in nodes",
            (averageu32(&paths_len)).to_string().as_str(),
        );
        file_info.add_info(
            "Longest path in bps",
            (*paths_bp_len.iter().max().unwrap()).to_string().as_str(),
        );
        file_info.add_info(
            "Shortest path in bps",
            (*paths_bp_len.iter().min().unwrap()).to_string().as_str(),
        );
        file_info.add_info(
            "Average path length in bps",
            (averageu32(&paths_bp_len)).to_string().as_str(),
        );
        file_info
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
    ) -> anyhow::Result<(
        ItemTable,
        Vec<String>,
        Vec<usize>,
        Vec<String>,
        HashMap<PathSegment, (u32, u32)>,
    )> {
        log::info!("parsing path + walk sequences");
        let (item_table, exclude_table, subset_covered_bps, paths_len, _collected_paths) =
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

        let path_names: Vec<String> = self
            .graph_storage
            .path_segments
            .iter()
            .map(|x| x.to_string())
            .collect();
        let group_names: Vec<String> = self
            .graph_storage
            .path_segments
            .iter()
            .map(|x| graph_mask.groups[x].clone())
            .collect_vec();
        let mgroup_names: HashSet<&str> = graph_mask
            .get_path_order(&self.graph_storage.path_segments)
            .into_iter()
            .map(|x| x.1)
            .collect();
        let (item_table, groups) = if group_names == path_names
            && group_names.len() == mgroup_names.len()
            && exclude_table.is_none()
        {
            (item_table, group_names)
        } else {
            Self::collapse_item_table(item_table, group_names, mgroup_names, &exclude_table)
        };
        let (feature_lengths, feature_names) =
            Self::get_feature_lengths(graph_storage, &exclude_table, &subset_covered_bps, count);

        Ok((
            item_table,
            groups,
            feature_lengths,
            feature_names,
            paths_len,
        ))
    }

    fn get_feature_lengths(
        graph_storage: &GraphStorage,
        exclude_table: &Option<ActiveTable>,
        subset_covered_bps: &Option<IntervalContainer>,
        count_type: CountType,
    ) -> (Vec<usize>, Vec<String>) {
        let mut feature_lengths = match count_type {
            CountType::Node => vec![1; graph_storage.node_count],
            CountType::Edge => vec![1; graph_storage.edge_count],
            CountType::Bp => graph_storage
                .node_lens
                .iter()
                .skip(1)
                .map(|x| *x as usize)
                .collect_vec(),
        };
        let mut feature_names = match count_type {
            CountType::Node | CountType::Bp => {
                let rev: HashMap<&ItemId, &str> = graph_storage
                    .node2id
                    .iter()
                    .map(|(k, v)| (v, std::str::from_utf8(k).unwrap()))
                    .collect();
                let names: Vec<String> = (1..feature_lengths.len() + 1)
                    .map(|x| rev[&ItemId(x as u64)].to_string())
                    .collect();
                names
            }
            CountType::Edge => {
                let rev: HashMap<&ItemId, String> = graph_storage
                    .edge2id
                    .as_ref()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (v, k.to_string()))
                    .collect();
                let names: Vec<String> = (1..feature_lengths.len() + 1)
                    .map(|x| rev[&ItemId(x as u64)].to_string())
                    .collect();
                names
            }
        };

        if let Some(e) = exclude_table.as_ref() {
            // Shorten length
            for annotated_feature in e.get_annotation_keys() {
                let idx = annotated_feature.0 as usize - 1;
                for (start, end) in e.get_active_intervals(annotated_feature, 0) {
                    let len = end - start;
                    feature_lengths[idx] -= len;
                }
            }
            // Iterator over excluded features (this has a zero element)
            // TODO can we just remove the zero element in the exclude_table
            // definition?
            let mut iter = e.items.iter().skip(1);
            feature_lengths.retain(|_| !iter.next().unwrap());

            // Do the same for the feature names
            let mut iter = e.items.iter().skip(1);
            feature_names.retain(|_| !iter.next().unwrap());
        }

        if let Some(s) = subset_covered_bps.as_ref() {
            if count_type == CountType::Bp {
                for key in s.keys() {
                    let feature_idx = key.0 as usize - 1;
                    let subset_length: usize = s
                        .get(key)
                        .unwrap()
                        .iter()
                        .map(|(start, end)| end - start)
                        .sum();
                    if subset_length < feature_lengths[feature_idx] {
                        feature_lengths[feature_idx] = subset_length;
                    }
                }
            }
        }

        (feature_lengths, feature_names)
    }

    fn collapse_item_table(
        item_table: ItemTable,
        groups: Vec<String>,
        masked_groups: HashSet<&str>,
        exclude_table: &Option<ActiveTable>,
    ) -> (ItemTable, Vec<String>) {
        let mut items: HashMap<&str, Vec<ItemIdSize>> = HashMap::new();
        // Translation table to compress id space (in case features are excluded)
        let translation_of_ids: Option<HashMap<usize, usize>> = exclude_table.as_ref().map(|e| {
            e.items
                .iter()
                .enumerate()
                .filter_map(|(idx, active)| if *active { None } else { Some(idx) })
                .enumerate()
                .map(|(new_idx, old_idx)| (old_idx, new_idx))
                .collect()
        });
        log::info!("Translation table: {:?}", translation_of_ids);
        for (idx, group) in groups.iter().enumerate() {
            let (start, end) = (
                item_table.id_prefsum[idx] as usize,
                item_table.id_prefsum[idx + 1] as usize,
            );
            items
                .entry(group)
                .or_insert(Vec::new())
                .extend(&item_table.items[start..end]);
        }
        let mut new_items = Vec::new();
        let mut new_idprefsum = vec![0];
        let mut keys = Vec::new();
        for (key, mut value) in items.into_iter() {
            // If path is excluded, skip it
            if !masked_groups.contains(key) {
                continue;
            }
            if let Some(e) = exclude_table.as_ref() {
                // If the feature should be removed, remove it
                value.retain(|x| !e.is_active(&ItemId(*x)));
                // Compress the id space (excluded features should take no index)
                value.iter_mut().for_each(|x| {
                    let new_idx = translation_of_ids.as_ref().unwrap()[&(*x as usize)];
                    *x = new_idx as u64;
                });
            }
            let len = value.len();
            new_items.append(&mut value);
            new_idprefsum.push(new_idprefsum.last().unwrap() + len as u64);
            keys.push(key.to_string());
        }
        (
            ItemTable {
                items: new_items,
                id_prefsum: new_idprefsum,
            },
            keys,
        )
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

        let mut fully_included = vec![false; graph_storage.node_count + 1];

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
                                &mut fully_included,
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

fn skip_path(item_table: &mut ItemTable, num_path: &mut usize, buf: &mut Vec<u8>) {
    item_table.id_prefsum[*num_path + 1] += item_table.id_prefsum[*num_path];
    *num_path += 1;
    buf.clear();
}

fn connected_components(edge2id: &HashMap<Edge, ItemId>, nodes: &Vec<ItemId>) -> Vec<u32> {
    let mut component_lengths = Vec::new();
    let mut visited: HashSet<ItemId> = HashSet::new();
    let edges: HashMap<ItemId, Vec<ItemId>> = edge2id
        .keys()
        .map(|x| (x.0, x.2))
        .chain(edge2id.keys().map(|x| (x.2, x.0)))
        .fold(HashMap::new(), |mut acc, (k, v)| {
            acc.entry(k).and_modify(|x| x.push(v)).or_insert(vec![v]);
            acc
        });
    for node in nodes {
        if !visited.contains(node) {
            component_lengths.push(dfs(&edges, *node, &mut visited));
        }
    }
    component_lengths
}

fn dfs(edges: &HashMap<ItemId, Vec<ItemId>>, node: ItemId, visited: &mut HashSet<ItemId>) -> u32 {
    let mut s = Vec::new();
    let mut length = 0;
    s.push(node);
    while let Some(v) = s.pop() {
        if visited.contains(&v) {
            continue;
        }
        visited.insert(v);
        length += 1;
        if !edges.contains_key(&v) {
            continue;
        }
        for neigh in &edges[&v] {
            if !visited.contains(neigh) {
                s.push(*neigh);
            }
        }
    }
    length
}
