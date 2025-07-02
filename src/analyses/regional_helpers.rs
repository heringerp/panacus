use std::collections::{HashMap, HashSet, VecDeque};

use crate::graph_broker::ItemId;

pub fn get_close_nodes(
    ref_nodes: &Vec<ItemId>,
    neighbors: &HashMap<ItemId, Vec<ItemId>>,
) -> HashMap<ItemId, Vec<ItemId>> {
    let mut queue: VecDeque<ItemId> = VecDeque::new();
    let mut closest_ref: HashMap<ItemId, ItemId> = HashMap::new();
    let ref_node_set: HashSet<ItemId> = ref_nodes.iter().copied().collect();
    for ref_node in ref_nodes {
        closest_ref.insert(*ref_node, *ref_node);
        for neighbor in neighbors.get(ref_node).unwrap_or(&Vec::new()) {
            if !ref_node_set.contains(neighbor) && !closest_ref.contains_key(neighbor) {
                closest_ref.insert(*neighbor, *ref_node);
                queue.push_front(*neighbor);
            }
        }
    }
    while let Some(current_node) = queue.pop_back() {
        for neighbor in neighbors.get(&current_node).unwrap_or(&Vec::new()) {
            if !closest_ref.contains_key(neighbor) {
                closest_ref.insert(*neighbor, current_node);
                queue.push_front(*neighbor);
            }
        }
    }
    let mut flipped: HashMap<ItemId, Vec<ItemId>> = HashMap::new();
    for (entry_node, ref_node) in closest_ref {
        flipped.entry(ref_node).or_default().push(entry_node);
    }
    flipped
}

pub fn get_ref_length(ref_nodes: &Vec<ItemId>, node_lens: &Vec<u32>) -> u32 {
    let length = ref_nodes
        .iter()
        .map(|node| node_lens[node.0 as usize])
        .sum();
    length
}

pub fn get_windows(
    ref_nodes: &Vec<ItemId>,
    node_lens: &Vec<u32>,
    number_of_windows: usize,
) -> (Vec<Vec<(ItemId, usize)>>, usize) {
    let ref_length = get_ref_length(ref_nodes, node_lens) as usize;
    let window_size = ref_length as usize / number_of_windows;
    let window_size = if window_size * number_of_windows == ref_length as usize {
        window_size
    } else {
        window_size + 1 // Add 1 to make last window uncomplete/total number stay at number_of_windows
    };
    let mut windows: Vec<Vec<(ItemId, usize)>> = vec![Vec::new(); number_of_windows];
    let mut bp_counter = 0;
    let mut current_window_index = 0;
    let mut current_node_index = 0;
    let mut already_used_bps_of_node = 0;
    while bp_counter < ref_length {
        let current_window_length: usize =
            windows[current_window_index].iter().map(|(_, l)| l).sum();
        let remaining_bps_in_window = window_size - current_window_length as usize;
        let current_node_length =
            node_lens[ref_nodes[current_node_index].0 as usize] as usize - already_used_bps_of_node;
        if current_node_length <= remaining_bps_in_window {
            windows[current_window_index]
                .push((ref_nodes[current_node_index], current_node_length as usize));
            bp_counter += current_node_length;
            current_node_index += 1;
            already_used_bps_of_node = 0;
            if current_node_length == remaining_bps_in_window {
                current_window_index += 1;
            }
        } else {
            let cut_node_length = remaining_bps_in_window;
            windows[current_window_index]
                .push((ref_nodes[current_node_index], cut_node_length as usize));
            bp_counter += cut_node_length;
            current_window_index += 1;
            already_used_bps_of_node += cut_node_length;
        };
    }
    (windows, window_size)
}
