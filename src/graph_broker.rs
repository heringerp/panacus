use core::panic;
use std::collections::VecDeque;
use std::io::BufRead;
use std::iter::zip;
use std::mem::take;
use std::time::Instant;
use std::{
    collections::{HashMap, HashSet},
    io::{BufWriter, Error, Write},
    str,
};

use abacus::{AbacusByTotal, GraphMask};
use graph::GraphStorage;
use util::{
    parse_path_identifier, parse_path_seq_to_item_vec, parse_path_seq_update_tables_multiple,
    parse_walk_identifier, parse_walk_seq_to_item_vec, parse_walk_seq_update_tables_multiple,
    update_tables_edgecount, update_tables_multiple,
};

use crate::util::{intersects, is_contained, ActiveTable, IntervalContainer, ItemTable, Threshold};
use crate::{
    analyses::InputRequirement as Req, analysis_parameter::Grouping,
    io::bufreader_from_compressed_gfa, util::CountType,
};

mod abacus;
mod coverage_matrix;
mod graph;
mod hist;
mod util;

pub use abacus::GraphMaskParameters;
pub use coverage_matrix::CoverageMatrix;
pub use graph::Edge;
pub use graph::ItemId;
pub use graph::Orientation;
pub use graph::PathSegment;
pub use hist::Hist;
pub use hist::ThresholdContainer;

#[derive(Debug, Clone, Default)]
pub struct GraphState {
    pub graph: String,
    pub name: Option<String>,
    pub subset: String,
    pub exclude: String,
    pub grouping: Option<Grouping>,
}

#[derive(Debug, Clone)]
pub struct GraphBroker {
    state: Option<GraphState>,
    // GraphStorage
    graph_aux: Option<GraphStorage>,
    name: String,

    // AbabcusAuxilliary
    abacus_aux_params: GraphMaskParameters,
    graph_mask: Option<GraphMask>,

    total_abaci: Option<HashMap<CountType, AbacusByTotal>>,
    group_abacus: Option<CoverageMatrix>,
    hists: Option<HashMap<CountType, Hist>>,
    csc_abacus: bool,
    paths: HashMap<PathSegment, Vec<(ItemId, Orientation)>>,

    path_lens: Option<HashMap<PathSegment, (u32, u32)>>,
    gfa_file: String,
    _nice: bool,
    input_requirements: HashSet<Req>,
    count_type: CountType,
}

impl GraphBroker {
    pub fn new() -> Self {
        GraphBroker {
            state: None,
            name: "".to_string(),
            graph_aux: None,
            abacus_aux_params: GraphMaskParameters::default(),
            graph_mask: None,
            total_abaci: None,
            group_abacus: None,
            hists: None,
            _nice: false,
            path_lens: None,
            gfa_file: String::new(),
            input_requirements: HashSet::new(),
            count_type: CountType::All,
            csc_abacus: false,
            paths: HashMap::new(),
        }
    }

    // TODO: fix situation instead of calculating the third value unnecessary
    fn contains_at_least_two(input_requirements: &HashSet<Req>) -> bool {
        if input_requirements.contains(&Req::Node) && input_requirements.contains(&Req::Edge) {
            return true;
        } else if input_requirements.contains(&Req::Edge) && input_requirements.contains(&Req::Bp) {
            return true;
        } else if input_requirements.contains(&Req::Bp) && input_requirements.contains(&Req::Node) {
            return true;
        } else {
            return false;
        }
    }

    pub fn change_graph_state(
        &mut self,
        state: GraphState,
        input_requirements: &HashSet<Req>,
        nice: bool,
    ) -> Result<(), Error> {
        if self.state.is_some() {
            let prev_state = std::mem::take(&mut self.state).unwrap();
            if prev_state.graph != state.graph {
                *self = Self::from_gfa(input_requirements, nice);
            }
            if prev_state.subset != state.subset {
                self.include_coords(&state.subset);
            }
            if prev_state.exclude != state.exclude {
                self.exclude_coords(&state.exclude);
            }
            if prev_state.grouping != state.grouping {
                self.with_group(&state.grouping);
            }
            if let Some(name) = &state.name {
                self.name = name.to_owned();
            } else {
                self.name = self.get_default_run_name(&state);
            }
            self.finish()?;
        } else {
            *self = Self::from_gfa(input_requirements, nice);
            if !state.subset.is_empty() {
                self.include_coords(&state.subset);
            }
            if !state.exclude.is_empty() {
                self.exclude_coords(&state.exclude);
            }
            if state.grouping.is_some() {
                self.with_group(&state.grouping);
            }
            if let Some(name) = &state.name {
                self.name = name.to_owned();
            } else {
                self.name = self.get_default_run_name(&state);
            }
            self.finish()?;
        }
        self.state = Some(state);
        Ok(())
    }

    pub fn change_order(&mut self, order: Option<&str>) -> Result<(), Error> {
        self.with_order(order);
        self.finish()
    }

    pub fn get_node_names(&self) -> HashMap<ItemId, Vec<u8>> {
        let tuples = self
            .graph_aux
            .as_ref()
            .expect("GraphBroker contains graph")
            .get_node_tuples();
        tuples.into_iter().map(|(name, id)| (id, name)).collect()
    }

    pub fn get_node_coverages(&self) -> HashMap<ItemId, u64> {
        self.group_abacus
            .as_ref()
            .expect("Contains group abacus")
            .get_node_coverages()
    }

    fn from_gfa(input_requirements: &HashSet<Req>, nice: bool) -> Self {
        let count_type = if Self::contains_at_least_two(input_requirements) {
            CountType::All
        } else if input_requirements.contains(&Req::Node) {
            CountType::Node
        } else if input_requirements.contains(&Req::Bp) {
            CountType::Bp
        } else if input_requirements.contains(&Req::Edge) {
            CountType::Edge
        } else {
            CountType::Node
        };
        let gfa_file = input_requirements
            .iter()
            .find(|v| matches!(v, Req::Graph(_)))
            .expect("Requirements contain gfa file");
        let gfa_file = match gfa_file {
            Req::Graph(gfa_file) => gfa_file,
            _ => panic!("Requirements really need to contain gfa file"),
        };
        let graph_aux = Some(GraphStorage::from_gfa(gfa_file, nice, count_type));
        GraphBroker {
            state: None,
            name: "".to_string(),
            graph_aux,
            abacus_aux_params: GraphMaskParameters::default(),
            graph_mask: None,
            total_abaci: None,
            group_abacus: None,
            hists: None,
            path_lens: None,
            gfa_file: gfa_file.to_owned(),
            _nice: nice,
            input_requirements: input_requirements.clone(),
            count_type,
            csc_abacus: false,
            paths: HashMap::new(),
        }
    }

    fn with_group(&mut self, grouping: &Option<Grouping>) {
        if let Some(grouping) = grouping {
            match grouping {
                Grouping::Sample => self.with_sample_group(),
                Grouping::Haplotype => self.with_haplo_group(),
                Grouping::Custom(file_name) => self.with_custom_group(file_name),
            };
        }
    }

    fn with_custom_group(&mut self, file_name: &str) {
        self.abacus_aux_params.groupby = file_name.to_owned();
    }

    fn with_haplo_group(&mut self) {
        self.abacus_aux_params.groupby_haplotype = true;
    }

    fn with_sample_group(&mut self) {
        self.abacus_aux_params.groupby_sample = true;
    }

    fn include_coords(&mut self, file_name: &str) {
        self.abacus_aux_params.positive_list = file_name.to_owned();
    }

    fn exclude_coords(&mut self, exclude: &str) {
        self.abacus_aux_params.negative_list = exclude.to_owned();
    }

    fn with_order(&mut self, file_name: Option<&str>) {
        self.abacus_aux_params.order = file_name.map(str::to_owned);
    }

    pub fn with_csc_abacus(mut self) -> Self {
        self.csc_abacus = true;
        self
    }

    // Second parse of graph (this time depending on subset, etc.)
    // focussed on graph paths
    pub fn finish(&mut self) -> Result<(), Error> {
        self.set_graph_mask()?;
        self.set_abaci_by_total();
        if self.input_requirements.contains(&Req::Hist) {
            self.set_hists();
        }
        let mut has_already_used_abacus = false;
        for req in self.input_requirements.clone() {
            match req {
                Req::AbacusByGroup(count) => {
                    if has_already_used_abacus {
                        panic!("Panacus is currently not able to have multiple Abaci By Group for different countables. Please run panacus either multiple times or wait for the planned pipelining feature");
                    }
                    self.set_abacus_by_group(count)?;
                    has_already_used_abacus = true;
                }
                _ => continue,
            }
        }
        Ok(())
    }

    pub fn get_run_name(&self) -> String {
        self.name.to_owned()
    }

    pub fn get_run_id(&self) -> String {
        self.name
            .to_lowercase()
            .replace([' ', '_', '#', '/', '"'], "-")
            .to_owned()
    }

    fn get_default_run_name(&self, state: &GraphState) -> String {
        if state.grouping.is_some() {
            format!(
                "{}-{}-{}",
                state.graph,
                state.subset,
                state.grouping.as_ref().unwrap()
            )
        } else {
            format!("{}-{}", state.graph, state.subset)
        }
    }

    pub fn get_degree(&self) -> &Vec<u32> {
        Self::check_and_error(self.graph_aux.as_ref().unwrap().degree.as_ref(), "degree");
        self.graph_aux.as_ref().unwrap().degree.as_ref().unwrap()
    }

    pub fn get_node_lens(&self) -> &Vec<u32> {
        &self.graph_aux.as_ref().unwrap().node_lens
    }

    pub fn get_edges(&self) -> &HashMap<Edge, ItemId> {
        Self::check_and_error(self.graph_aux.as_ref().unwrap().edge2id.as_ref(), "edge2id");
        self.graph_aux.as_ref().unwrap().edge2id.as_ref().unwrap()
    }

    pub fn get_nodes(&self) -> Vec<ItemId> {
        self.graph_aux.as_ref().unwrap().get_nodes()
    }

    pub fn get_node_count(&self) -> usize {
        self.graph_aux.as_ref().unwrap().node_count
    }

    pub fn get_edge_count(&self) -> usize {
        self.graph_aux.as_ref().unwrap().edge_count
    }

    pub fn get_group_count(&self) -> usize {
        Self::check_and_error(self.graph_mask.as_ref(), "abacus_aux -> group_count");
        self.graph_mask.as_ref().unwrap().count_groups()
    }

    pub fn get_fname(&self) -> String {
        self.gfa_file.to_string()
    }

    pub fn get_groups(&self) -> &HashMap<PathSegment, String> {
        Self::check_and_error(self.graph_mask.as_ref(), "abacus_aux -> groups");
        &self.graph_mask.as_ref().unwrap().groups
    }

    pub fn get_path_lens(&self) -> &HashMap<PathSegment, (u32, u32)> {
        Self::check_and_error(self.path_lens.as_ref(), "path_lens");
        self.path_lens.as_ref().unwrap()
    }

    pub fn get_hists(&self) -> &HashMap<CountType, Hist> {
        Self::check_and_error(self.hists.as_ref(), "hists");
        self.hists.as_ref().unwrap()
    }

    pub fn get_abacus_by_group(&self) -> &CoverageMatrix {
        Self::check_and_error(self.group_abacus.as_ref(), "abacus_by_group");
        self.group_abacus.as_ref().unwrap()
    }

    pub fn get_abacus_by_total(&self, count: CountType) -> &AbacusByTotal {
        Self::check_and_error(self.total_abaci.as_ref(), "abacus_by_group");
        &self.total_abaci.as_ref().unwrap()[&count]
    }

    pub fn get_growth_for_subset(
        &self,
        count: CountType,
        indices: &Vec<usize>,
        uncovered_bps: Option<HashMap<u64, usize>>,
        coverage: usize,
    ) -> Vec<f64> {
        let abacus = self.get_abacus_by_total(count);
        let hist =
            Hist::from_abacus_for_window(abacus, self.graph_aux.as_ref(), indices, uncovered_bps);
        let cov = Threshold::Absolute(coverage);
        let growth = hist.calc_growth_union(&cov);
        growth
    }

    pub fn get_hist_for_subset(
        &self,
        count: CountType,
        indices: &Vec<usize>,
        uncovered_bps: Option<HashMap<u64, usize>>,
    ) -> Vec<f64> {
        let abacus = self.get_abacus_by_total(count);
        let hist =
            Hist::from_abacus_for_window(abacus, self.graph_aux.as_ref(), indices, uncovered_bps);
        hist.coverage.into_iter().map(|c| c as f64).collect()
    }

    pub fn get_hist_for_path_subset(
        &self,
        count: CountType,
        path_indices: &Vec<usize>,
        uncovered_bps: Option<HashMap<u64, usize>>,
    ) -> Vec<f64> {
        Vec::new()
    }

    pub fn get_path(&self, path_seg: &PathSegment) -> &Vec<(ItemId, Orientation)> {
        &self.paths[path_seg]
    }

    pub fn get_all_matchings_paths(
        &self,
        path_seg: &PathSegment,
    ) -> HashMap<PathSegment, &Vec<(ItemId, Orientation)>> {
        let mut result = HashMap::new();
        for path in self.paths.keys() {
            if path.is_part_of(path_seg) {
                result.insert(path.clone(), &self.paths[path]);
            }
        }
        result
    }

    pub fn write_abacus_by_group<W: Write>(
        &self,
        total: bool,
        out: &mut BufWriter<W>,
    ) -> Result<(), Error> {
        Self::check_and_error(self.group_abacus.as_ref(), "abacus_by_group");
        self.group_abacus
            .as_ref()
            .unwrap()
            .to_tsv(total, out, self.graph_aux.as_ref().unwrap())
    }

    fn set_graph_mask(&mut self) -> Result<(), Error> {
        self.graph_mask = Some(GraphMask::from_datamgr(
            &self.abacus_aux_params,
            self.graph_aux.as_ref().unwrap(),
        )?);
        Ok(())
    }

    fn set_hists(&mut self) {
        let mut hists = HashMap::new();
        for (k, v) in self.total_abaci.as_ref().unwrap() {
            hists.insert(
                *k,
                Hist::from_abacus(v, Some(self.graph_aux.as_ref().unwrap())),
            );
        }
        self.hists = Some(hists);
    }

    fn check_and_error<T>(value: Option<T>, type_of_value: &str) {
        if value.is_none() {
            let msg = format!(
                "Cannot give value of {}, since it was not requested",
                type_of_value
            );
            log::error!("{}", &msg);
        }
    }

    fn set_abacus_by_group(&mut self, count: CountType) -> Result<(), Error> {
        // let mut abaci_by_group = HashMap::new();
        let mut data = bufreader_from_compressed_gfa(&self.gfa_file);
        let abacus = CoverageMatrix::from_gfa(
            &mut data,
            self.graph_mask.as_ref().unwrap(),
            self.graph_aux.as_ref().unwrap(),
            count,
            true,
        )?;
        // abaci_by_group.insert(self.count_type, abacus);
        self.group_abacus = Some(abacus);
        Ok(())
    }

    fn get_paths_to_collect(&self) -> Vec<PathSegment> {
        self.input_requirements
            .iter()
            .filter_map(|r| match r {
                Req::Path(p) => Some(p.clone()),
                _ => None,
            })
            .collect()
    }

    fn set_abaci_by_total(&mut self) {
        let (count_types_not_edge, shall_calculate_edge) = self.get_abacus_count_type();
        log::info!(
            "calculating abaci for count_types: {:?}, and edge: {}",
            count_types_not_edge,
            shall_calculate_edge
        );
        let mut abaci = if !count_types_not_edge.is_empty() {
            let paths_to_collect = self.get_paths_to_collect();
            let (abaci, path_lens, collected_paths) =
                self.set_abacus_from_gfa(&count_types_not_edge, &paths_to_collect);
            let abaci: HashMap<CountType, AbacusByTotal> =
                zip(count_types_not_edge, abaci).collect();
            if self.input_requirements.contains(&Req::PathLens) {
                self.path_lens = Some(path_lens);
            }
            self.paths = collected_paths;
            abaci
        } else {
            HashMap::new()
        };
        if shall_calculate_edge {
            let (mut edge_abacus, _, _) =
                self.set_abacus_from_gfa(&vec![CountType::Edge], &Vec::new());
            abaci.insert(CountType::Edge, edge_abacus.pop().unwrap());
        }
        self.total_abaci = Some(abaci);
    }

    // Nodes and bps behave similarly, while edges have to be calculated separately.
    fn get_abacus_count_type(&mut self) -> (Vec<CountType>, bool) {
        let count_types_not_edge = if self.count_type == CountType::All {
            vec![CountType::Node, CountType::Bp]
        } else if self.count_type != CountType::Edge {
            vec![self.count_type.clone()]
        } else {
            Vec::new()
        };
        let shall_calculate_edge =
            self.count_type == CountType::All || self.count_type == CountType::Edge;
        (count_types_not_edge, shall_calculate_edge)
    }

    pub fn get_hist_from_abacus(&self, abacus: &AbacusByTotal) -> Hist {
        Hist::from_abacus(abacus, self.graph_aux.as_ref())
    }

    pub fn set_abacus_from_gfa(
        &self,
        count_types: &Vec<CountType>,
        paths_to_collect: &Vec<PathSegment>,
    ) -> (
        Vec<AbacusByTotal>,
        HashMap<PathSegment, (u32, u32)>,
        HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    ) {
        let (item_tables, exclude_tables, mut subset_covered_bps, path_lens, collected_paths) =
            self.parse_paths_walks(count_types, paths_to_collect);
        let mut item_tables = VecDeque::from(item_tables);
        let mut exclude_tables = VecDeque::from(exclude_tables);
        let mut subset_covered_bps: VecDeque<_> = count_types
            .iter()
            .map(|count| match count {
                &CountType::Bp if subset_covered_bps.is_some() => take(&mut subset_covered_bps),
                _ => None,
            })
            .collect();
        let abaci = count_types
            .iter()
            .map(|count| {
                AbacusByTotal::item_table_to_abacus(
                    self.graph_mask
                        .as_ref()
                        .expect("Graph mask should exist at this point"),
                    self.graph_aux
                        .as_ref()
                        .expect("Graph storage should exist at this point"),
                    *count,
                    item_tables.pop_front().unwrap(),
                    exclude_tables.pop_front().unwrap(),
                    subset_covered_bps.pop_front().unwrap(),
                )
            })
            .collect();
        (abaci, path_lens, collected_paths)
    }

    pub fn parse_paths_walks(
        &self,
        count_types: &Vec<CountType>,
        paths_to_collect: &Vec<PathSegment>,
    ) -> (
        Vec<ItemTable>,
        Vec<Option<ActiveTable>>,
        Option<IntervalContainer>,
        HashMap<PathSegment, (u32, u32)>,
        HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    ) {
        log::info!("parsing path + walk sequences");
        log::info!("collecting: {:?}", paths_to_collect);
        let mut data = bufreader_from_compressed_gfa(&self.gfa_file);
        let graph_storage = self
            .graph_aux
            .as_ref()
            .expect("Graph should have set graph storage");
        let graph_mask = self
            .graph_mask
            .as_ref()
            .expect("Graph should have graph mask set");
        let mut item_tables =
            vec![ItemTable::new(graph_storage.path_segments.len()); count_types.len()];

        let (mut subset_covered_bps, mut exclude_tables, include_map, exclude_map) =
            graph_mask.load_optional_subsetting_multiple(graph_storage, count_types);

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
                    skip_path(&mut item_tables, &mut num_path, &mut buf);
                    continue;
                }

                let indices = get_count_indices(count_types);
                indices.into_iter().for_each(|(count, is)| {
                if count != CountType::Edge
                    && (graph_mask.include_coords.is_none()
                        || is_contained(include_coords, &(start, end)))
                        && (graph_mask.exclude_coords.is_none()
                            || is_contained(exclude_coords, &(start, end)))
                {
                    log::debug!("path {} is fully contained within subset coordinates {:?} and is eligible for full parallel processing", path_seg, include_coords);
                    let mut none = None;
                    let ex: Vec<&mut Option<ActiveTable>> = if exclude_coords.is_empty() {
                        vec![&mut none]
                    } else {
                        exclude_tables.iter_mut().enumerate().filter(|(i, _)| is.contains(i)).map(|(_, e)| e).collect()
                    };
                    let (num_added_nodes, bp_len, path_sequence) = match buf[0] {
                        b'P' => parse_path_seq_update_tables_multiple(
                            buf_path_seg,
                            graph_storage,
                            &mut item_tables[is[0]],
                            ex,
                            num_path,
                            paths_to_collect.iter().any(|p| path_seg.is_part_of(p))
                        ),
                        b'W' => parse_walk_seq_update_tables_multiple(
                            buf_path_seg,
                            graph_storage,
                            &mut item_tables[is[0]],
                            ex,
                            num_path,
                            paths_to_collect.iter().any(|p| path_seg.is_part_of(p))
                        ),
                        _ => unreachable!(),
                    };
                    if let Some(path_sequence) = path_sequence {
                        log::info!("inserting {:?}", path_seg);
                        collected_paths.insert(path_seg.clone(), path_sequence);
                    }
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
                    let mut exclude_tables_red = exclude_tables.iter_mut().enumerate().filter(|(i, _)| is.contains(i)).map(|(_, e)| e).collect();
                    match count {
                        CountType::Node | CountType::Bp => {
                            let (node_len, bp_len) = update_tables_multiple(
                                &mut item_tables[is[0]],
                                &mut subset_covered_bps.as_mut(),
                                exclude_tables_red,
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
                            &mut item_tables[is[0]],
                            &mut exclude_tables_red[0].as_mut(),
                            num_path,
                            graph_storage,
                            sids,
                            include_coords,
                            exclude_coords,
                            start,
                        ),
                        CountType::All => unreachable!("inadmissable count type"),
                    };
                }
            });
                num_path += 1;
            }
            buf.clear();
        }
        let duration = timer.elapsed();
        log::info!(
            "func done; count: {:?}; time elapsed: {:?}",
            count_types,
            duration
        );

        // In case of node/bps, copy first table into the second
        if item_tables.len() == 2 {
            item_tables[1] = item_tables[0].clone();
        }
        (
            item_tables,
            exclude_tables,
            subset_covered_bps,
            paths_len,
            collected_paths,
        )
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

fn skip_path(item_tables: &mut Vec<ItemTable>, num_path: &mut usize, buf: &mut Vec<u8>) {
    for item_table in item_tables {
        item_table.id_prefsum[*num_path + 1] += item_table.id_prefsum[*num_path];
    }
    *num_path += 1;
    buf.clear();
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
