use std::{
    collections::HashMap,
    io::{BufWriter, Error, Write},
};

use itertools::Itertools;

use crate::{
    graph_broker::{
        abacus::quantify_uncovered_bps, util::parse_gfa_paths_walks, Edge, Orientation,
    },
    util::{ActiveTable, CountSize, CountType, GroupSize, ItemIdSize, ItemTable, Threshold, Wrap},
};

use super::{
    abacus::{AbacusByTotal, GraphMask},
    graph::GraphStorage,
    ItemId,
};

#[derive(Debug, Clone)]
pub struct CoverageMatrix {
    count: CountType,
    r: Vec<usize>,
    v: Option<Vec<CountSize>>,
    c: Vec<GroupSize>,
    uncovered_bps: HashMap<ItemIdSize, usize>,
    groups: Vec<String>,
}

impl CoverageMatrix {
    pub fn get_count(&self) -> CountType {
        self.count
    }

    pub fn get_groups(&self) -> &Vec<String> {
        &self.groups
    }

    // TODO remove this (need to change analysis)
    pub fn get_r(&self) -> &Vec<usize> {
        &self.r
    }

    // TODO remove this (need to change analysis)
    pub fn get_c(&self) -> &Vec<GroupSize> {
        &self.c
    }

    pub fn get_nodes_of_paths(&self, paths: &Vec<usize>) -> Vec<usize> {
        let mut indices = Vec::new();
        for (i, c) in self.c.iter().enumerate() {
            if paths.contains(&(*c as usize)) {
                indices.push(i);
            }
        }
        let mut nodes = Vec::new();
        'outer: for (i, (s, e)) in self.r.iter().tuple_windows().enumerate() {
            for index in indices.iter() {
                if index >= s && index < e {
                    nodes.push(i);
                    // As soon as node is added, continue, don't try to add that node again
                    continue 'outer;
                }
            }
        }
        nodes
    }

    pub fn to_abacus_all_nodes(&self, paths: &Vec<usize>) -> (AbacusByTotal, Vec<usize>) {
        let no_nodes = self.r.len() - 1;
        let nodes: Vec<usize> = (0..no_nodes).collect();
        self.to_abacus(&nodes, paths)
    }

    pub fn to_abacus(&self, nodes: &Vec<usize>, paths: &Vec<usize>) -> (AbacusByTotal, Vec<usize>) {
        let (coverage_list, non_zero_nodes): (Vec<u32>, Vec<Option<usize>>) = nodes
            .iter()
            .map(|node| {
                let mut coverage_of_node = 0u32;
                let interval_start = self.r[*node];
                let interval_end = self.r[node + 1];
                for i in interval_start..interval_end {
                    if paths.contains(&(self.c[i] as usize)) {
                        coverage_of_node += 1;
                    }
                }
                (
                    coverage_of_node,
                    if coverage_of_node == 0 {
                        None
                    } else {
                        Some(*node)
                    },
                )
            })
            .unzip();
        let groups: Vec<String> = paths.iter().map(|p| self.groups[*p].clone()).collect();
        let non_zero_nodes: Vec<usize> = non_zero_nodes.into_iter().filter_map(|x| x).collect();
        (
            AbacusByTotal {
                count: self.count,
                countable: coverage_list,
                uncovered_bps: Some(self.uncovered_bps.clone()),
                groups,
            },
            non_zero_nodes,
        )
    }

    pub fn from_gfa<R: std::io::Read>(
        data: &mut std::io::BufReader<R>,
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        count: CountType,
        report_values: bool,
    ) -> Result<Self, Error> {
        log::info!("parsing path + walk sequences");
        let (item_table, exclude_table, subset_covered_bps, _paths_len) =
            parse_gfa_paths_walks(data, graph_mask, graph_storage, &count);

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

        let r = CoverageMatrix::compute_row_storage_space(
            &item_table,
            &exclude_table,
            &path_order,
            graph_storage.number_of_items(&count),
        );
        let (v, c) =
            CoverageMatrix::compute_column_values(&item_table, &path_order, &r, report_values);
        log::info!(
            "abacus has {} path groups and {} countables",
            groups.len(),
            r.len()
        );
        let uncovered_bps =
            quantify_uncovered_bps(&exclude_table, &subset_covered_bps, graph_storage);

        Ok(Self {
            count,
            r,
            v,
            c,
            uncovered_bps,
            groups,
        })
    }

    fn compute_row_storage_space(
        item_table: &ItemTable,
        exclude_table: &Option<ActiveTable>,
        path_order: &Vec<(ItemIdSize, GroupSize)>,
        n_items: usize,
    ) -> Vec<usize> {
        log::info!("computing space allocating storage for group-based coverage table:");
        let mut last: Vec<GroupSize> = vec![GroupSize::MAX; n_items + 1];
        let last_ptr = Wrap(&mut last);

        let mut r: Vec<usize> = vec![0; n_items + 2];
        let r_ptr = Wrap(&mut r);
        for (path_id, group_id) in path_order {
            let start = item_table.id_prefsum[*path_id as usize] as usize;
            let end = item_table.id_prefsum[*path_id as usize + 1] as usize;
            for j in start..end {
                let sid = item_table.items[j] as usize;
                if &last[sid] != group_id
                    && (exclude_table.is_none() || !exclude_table.as_ref().unwrap().items[sid])
                {
                    unsafe {
                        (&mut (*r_ptr.0))[sid] += 1;
                        (&mut (*last_ptr.0))[sid] = *group_id;
                    }
                }
            }
        }
        log::info!(" ++ assigning storage locations");
        let mut c = 0;
        // can this be simplified?
        for item in &mut r {
            let tmp = *item;
            *item = c;
            c += tmp;
        }
        log::info!(
            " ++ group-aware table has {} non-zero elements",
            r.last().unwrap()
        );
        r
    }

    fn compute_column_values(
        item_table: &ItemTable,
        path_order: &Vec<(ItemIdSize, GroupSize)>,
        r: &[usize],
        report_values: bool,
    ) -> (Option<Vec<CountSize>>, Vec<GroupSize>) {
        let n = { *r.last().unwrap() };
        log::info!("allocating storage for group-based coverage table..");
        let mut v = if report_values {
            vec![0; n]
        } else {
            // we produce a dummy
            vec![0; 1]
        };
        let mut c: Vec<GroupSize> = vec![GroupSize::MAX; n];
        log::info!("done");

        log::info!("computing group-based coverage..");
        let v_ptr = Wrap(&mut v);
        let c_ptr = Wrap(&mut c);

        // group id is monotone increasing from 0 to #groups
        for (path_id, group_id) in path_order {
            let path_id_u = *path_id as usize;
            let start = item_table.id_prefsum[path_id_u] as usize;
            let end = item_table.id_prefsum[path_id_u + 1] as usize;
            for j in start..end {
                let sid = item_table.items[j] as usize;
                let cv_start = r[sid];
                let mut cv_end = r[sid + 1];
                if cv_end != cv_start {
                    // look up storage location for node cur_sid: we use the last position
                    // of interval cv_start..cv_end, which is associated to coverage counts
                    // of the current node (sid), in the "c" array as pointer to the
                    // current column (group) / value (coverage) position. If the current group
                    // id does not match the one associated with the current position, we move
                    // on to the next. If cv_start + p == cv_end - 1, this means that we are
                    // currently writing the last element in that interval, and we need to make
                    // sure that we are no longer using it as pointer.
                    if cv_end - 1 > c.len() {
                        log::error!(
                            "oops, cv_end-1 is larger than the length of c for sid={}",
                            sid
                        );
                        cv_end = c.len() - 1;
                    }

                    let mut p = c[cv_end - 1] as usize;
                    unsafe {
                        // we  look at an untouched interval, so let's get the pointer game
                        // started...
                        if c[cv_end - 1] == GroupSize::MAX {
                            (&mut (*c_ptr.0))[cv_start] = *group_id;
                            // if it's just a single value in this interval, the pointer game
                            // ends before it started
                            if cv_start < cv_end - 1 {
                                (&mut (*c_ptr.0))[cv_end - 1] = 0;
                            }
                            if report_values {
                                (&mut (*v_ptr.0))[cv_start] += 1;
                            }
                        } else if cv_start + p < cv_end - 1 {
                            // if group id of current slot does not match current group id
                            // (remember group id's are strictly monotically increasing), then
                            // move on to the next slot
                            if c[cv_start + p] < *group_id {
                                // move on to the next slot
                                (&mut (*c_ptr.0))[cv_end - 1] += 1;
                                // update local pointer
                                p += 1;
                                (&mut (*c_ptr.0))[cv_start + p] = *group_id
                            }
                            if report_values {
                                (&mut (*v_ptr.0))[cv_start + p] += 1;
                            }
                        } else if report_values {
                            // make sure it points to the last element and not beyond
                            (&mut (*v_ptr.0))[cv_end - 1] += 1;
                        }
                    }
                }
            }
        }
        log::info!("done");
        (if report_values { Some(v) } else { None }, c)
    }

    // why &self and not self? we could destroy abacus at this point.
    pub fn calc_growth(
        &self,
        t_coverage: &Threshold,
        t_quorum: &Threshold,
        node_lens: &Vec<u32>,
    ) -> Vec<f64> {
        let mut res = vec![0.0; self.groups.len()];

        let c = usize::max(1, t_coverage.to_absolute(self.groups.len()));
        let q = f64::max(0.0, t_quorum.to_relative(self.groups.len()));

        let mut it = self.r.iter().tuple_windows().enumerate();
        // ignore first entry
        it.next();
        for (i, (&start, &end)) in it {
            if end - start >= c {
                let mut k = start;
                for j in self.c[start] as usize..self.groups.len() {
                    if k < end - 1 && self.c[k + 1] as usize <= j {
                        k += 1
                    }
                    if k - start + 1 >= ((self.c[k] as f64 + 1.0) * q).ceil() as usize {
                        // we never need to look into the actual value in self.v, because we
                        // know it must be non-zero, which is sufficient
                        match self.count {
                            CountType::Node | CountType::Edge => res[j] += 1.0,
                            CountType::Bp => {
                                let uncovered =
                                    self.uncovered_bps.get(&(i as ItemIdSize)).unwrap_or(&0);
                                let covered = node_lens[i] as usize;
                                if uncovered > &covered {
                                    log::error!("oops, #uncovered bps ({}) is larger than #coverd bps ({}) for node with sid {})", &uncovered, &covered, i);
                                } else {
                                    res[j] += (covered - uncovered) as f64
                                }
                            }
                            CountType::All => unreachable!("inadmissible count type"),
                        }
                    }
                }
            }
        }
        res
    }

    #[allow(dead_code)]
    pub fn write_rcv<W: Write>(&self, out: &mut BufWriter<W>) -> Result<(), Error> {
        write!(out, "{}", self.r[0])?;
        for x in self.r[1..].iter() {
            write!(out, "\t{}", x)?;
        }
        writeln!(out)?;
        write!(out, "{}", self.c[0])?;
        for x in self.c[1..].iter() {
            write!(out, "\t{}", x)?;
        }
        writeln!(out)?;
        if let Some(v) = &self.v {
            write!(out, "{}", v[0])?;
            for x in v[1..].iter() {
                write!(out, "\t{}", x)?;
            }
            writeln!(out)?;
        };
        Ok(())
    }

    pub fn get_node_coverages(&self) -> HashMap<ItemId, u64> {
        let mut result = HashMap::new();

        match self.count {
            CountType::Node => {
                let mut it = self.r.iter().tuple_windows().enumerate();
                // ignore first entry
                it.next();
                for (i, (&start, &end)) in it {
                    let id = ItemId(i as ItemIdSize);
                    let total = end - start;
                    result.insert(id, total as u64);
                }
            }
            CountType::Edge | CountType::Bp | CountType::All => {
                unreachable!("inadmissible count type")
            }
        };

        result
    }

    pub fn to_tsv<W: Write>(
        &self,
        total: bool,
        out: &mut BufWriter<W>,
        graph_storage: &GraphStorage,
    ) -> Result<(), Error> {
        // create mapping from numerical node ids to original node identifiers
        log::info!("reporting coverage table");
        let dummy = Vec::new();
        let mut id2node: Vec<&Vec<u8>> = vec![&dummy; graph_storage.node_count + 1];
        // TODO: fix issue
        let node_tuples = graph_storage.get_node_tuples();
        for (node, id) in node_tuples.iter() {
            id2node[id.0 as usize] = node;
        }

        match self.count {
            CountType::Node | CountType::Bp => {
                write!(out, "node")?;
                if total {
                    write!(out, "\ttotal")?;
                } else {
                    for group in self.groups.iter() {
                        write!(out, "\t{}", group)?;
                    }
                }
                writeln!(out)?;

                let mut it = self.r.iter().tuple_windows().enumerate();
                // ignore first entry
                it.next();
                for (i, (&start, &end)) in it {
                    let bp = if self.count == CountType::Bp {
                        graph_storage.node_lens[i] as usize
                            - *self.uncovered_bps.get(&(i as ItemIdSize)).unwrap_or(&0)
                    } else {
                        1
                    };
                    write!(out, "{}", std::str::from_utf8(id2node[i]).unwrap())?;
                    if total {
                        // we never need to look into the actual value in self.v, because we
                        // know it must be non-zero, which is sufficient
                        writeln!(out, "\t{}", end - start)?;
                    } else {
                        let mut k = start;
                        for j in 0 as GroupSize..self.groups.len() as GroupSize {
                            if k == end || j < self.c[k] {
                                write!(out, "\t0")?;
                            } else if j == self.c[k] {
                                match &self.v {
                                    None => write!(out, "\t{}", bp),
                                    Some(v) => write!(out, "\t{}", v[k] as usize * bp),
                                }?;
                                k += 1;
                            }
                        }
                        writeln!(out)?;
                    }
                }
            }
            CountType::Edge => {
                if let Some(edge2id) = &graph_storage.edge2id {
                    let dummy_edge = crate::graph_broker::Edge(
                        ItemId(0),
                        Orientation::default(),
                        ItemId(0),
                        Orientation::default(),
                    );
                    let mut id2edge: Vec<&Edge> = vec![&dummy_edge; graph_storage.edge_count + 1];
                    for (edge, id) in edge2id.iter() {
                        id2edge[id.0 as usize] = edge;
                    }

                    write!(out, "edge")?;
                    if total {
                        write!(out, "\ttotal")?;
                    } else {
                        for group in self.groups.iter() {
                            write!(out, "\t{}", group)?;
                        }
                    }
                    writeln!(out)?;

                    let mut it = self.r.iter().tuple_windows().enumerate();
                    // ignore first entry
                    it.next();
                    for (i, (&start, &end)) in it {
                        let edge = id2edge[i];
                        write!(
                            out,
                            "{}{}{}{}",
                            edge.1,
                            std::str::from_utf8(id2node[edge.0 .0 as usize]).unwrap(),
                            edge.3,
                            std::str::from_utf8(id2node[edge.2 .0 as usize]).unwrap(),
                        )?;
                        if total {
                            // we never need to look into the actual value in self.v, because we
                            // know it must be non-zero, which is sufficient
                            writeln!(out, "\t{}", end - start)?;
                        } else {
                            let mut k = start;
                            for j in 0 as GroupSize..self.groups.len() as GroupSize {
                                if k == end || j < self.c[k] {
                                    write!(out, "\t0")?;
                                } else if j == self.c[k] {
                                    match &self.v {
                                        None => write!(out, "\t1"),
                                        Some(v) => write!(out, "\t{}", v[j as usize]),
                                    }?;
                                    k += 1;
                                }
                            }
                            writeln!(out)?;
                        }
                    }
                }
            }
            CountType::All => unreachable!("inadmissible count type"),
        };

        Ok(())
    }
}

//pub enum Abacus<'a> {
//    Total(AbacusByTotal<'a>),
//    Group(AbacusByGroup<'a>),
//    Nil,
//}
