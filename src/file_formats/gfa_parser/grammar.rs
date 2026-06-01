use lazy_static::lazy_static;
use std::{collections::HashSet, io::BufRead};

use regex::Regex;

use crate::{
    file_formats::gfa_parser::{graph::GraphStorage, ItemId, Orientation},
    io::bufreader_from_compressed_gfa,
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

pub struct Grammar {
    rules: Vec<usize>,
    rule_texts: Vec<(ItemId, Orientation)>,
}

impl Default for Grammar {
    fn default() -> Self {
        Self {
            rules: vec![0],
            rule_texts: Vec::new(),
        }
    }
}

impl Grammar {
    pub fn insert(&mut self, rule: usize, rule_text: Vec<(ItemId, Orientation)>) {
        if rule != self.rules.len() - 1 {
            log::error!("Grammar rules should be inserted in order");
        }
        self.rules
            .push(self.rules.last().unwrap() + rule_text.len());
        self.rule_texts.extend(rule_text);
    }

    fn get_topological_sort(&self, node2rule_id: &[usize]) -> Vec<ItemId> {
        log::info!("Calculating Topological Sort of Grammar");
        let number_of_rules = self.rules.len() - 1;
        let mut l: Vec<ItemId> = Vec::with_capacity(number_of_rules);
        let mut marked = vec![false; number_of_rules];
        for (idx, &rule_id) in node2rule_id.iter().enumerate() {
            if rule_id == usize::MAX {
                continue;
            }
            if !marked[rule_id] {
                let item_id = ItemId(idx as u64);
                self.visit(rule_id, item_id, &mut marked, &mut l, node2rule_id);
            }
        }
        l.reverse();
        l
    }

    fn visit(
        &self,
        rule: usize,
        item_id: ItemId,
        marked: &mut [bool],
        l: &mut Vec<ItemId>,
        node2rule_id: &[usize],
    ) {
        if marked[rule] {
            return;
        }
        marked[rule] = true;

        let (start, end) = (self.rules[rule], self.rules[rule + 1]);
        for (child, _) in &self.rule_texts[start..end] {
            if node2rule_id[child.0 as usize] != usize::MAX {
                self.visit(
                    node2rule_id[child.0 as usize],
                    *child,
                    marked,
                    l,
                    node2rule_id,
                );
            }
        }
        l.push(item_id)
    }

    pub fn push_down_values(&self, values: &mut Vec<u32>, node2rule_id: &[usize]) {
        let sorting = self.get_topological_sort(node2rule_id);
        // check topological sorting:
        let mut seen: HashSet<ItemId> = HashSet::new();
        for rule in sorting.iter().rev() {
            seen.insert(*rule);
            let rule_id = node2rule_id[rule.0 as usize];
            let (start, end) = (self.rules[rule_id], self.rules[rule_id + 1]);
            for (child, _) in &self.rule_texts[start..end] {
                let child_id = child.0 as usize;
                if node2rule_id[child_id] != usize::MAX && !seen.contains(child) {
                    log::error!(
                        "This non-terminal {:?} was not seen before, topological sort is wrong!",
                        child
                    );
                    panic!();
                }
            }
        }
        for rule in sorting {
            let value = values[rule.0 as usize];
            let rule_id = node2rule_id[rule.0 as usize];
            values[rule.0 as usize] = 0;
            let (start, end) = (self.rules[rule_id], self.rules[rule_id + 1]);
            for (child, _) in &self.rule_texts[start..end] {
                let child_id = child.0 as usize;
                values[child_id] += value;
            }
        }
    }

    pub fn get(&self, rule: usize) -> &[(ItemId, Orientation)] {
        let (start, end) = (self.rules[rule], self.rules[rule + 1]);
        &self.rule_texts[start..end]
    }

    pub fn parse_gfa(&mut self, filename: &str, graph_storage: &GraphStorage) {
        let buf = bufreader_from_compressed_gfa(filename);
        self.parse_buffer(buf, graph_storage);
    }

    pub fn is_empty(&self) -> bool {
        self.rule_texts.is_empty()
    }

    pub fn len(&self) -> usize {
        self.rules.len() - 1
    }

    fn parse_buffer(&mut self, buf: impl BufRead, graph_storage: &GraphStorage) {
        let lines = buf.lines().map(|l| l.expect("Failed to read line"));

        for line in lines {
            if line.starts_with("Q\t") {
                let mut fields = line.split("\t").skip(1);
                let name: Vec<u8> = fields
                    .next()
                    .expect("Q line has segment name")
                    .as_bytes()
                    .iter()
                    .copied()
                    .collect();
                let node_id = graph_storage.get_node_id(&name).unwrap();
                let rule_id = graph_storage.node2rule_id[node_id.0 as usize];
                let content = fields.next().expect("Q line has content");
                let nodes: Vec<(ItemId, Orientation)> = RE_WALK
                    .captures_iter(content)
                    .map(|m| {
                        let orientation = Orientation::from_lg(m[1].bytes().nth(0).unwrap() as u8);
                        let node: Vec<u8> = m[2].bytes().collect();
                        let node = graph_storage.get_node_id(&node).unwrap();
                        (node, orientation)
                    })
                    .collect();
                self.insert(rule_id, nodes);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let grammar = Grammar::default();
        assert_eq!(grammar.rules.len(), 1);
        assert_eq!(grammar.rule_texts.len(), 0);
    }

    #[test]
    fn test_insert_rule() {
        let mut grammar = Grammar::default();
        grammar.insert(
            0,
            vec![
                (ItemId(1), Orientation::Forward),
                (ItemId(2), Orientation::Backward),
            ],
        );
        assert_eq!(grammar.rules.len(), 2);
        assert_eq!(grammar.rule_texts.len(), 2);
    }

    #[test]
    fn test_lookup_rule() {
        let mut grammar = Grammar::default();
        let text = vec![
            (ItemId(1), Orientation::Forward),
            (ItemId(2), Orientation::Backward),
        ];
        grammar.insert(0, text.clone());
        let result = grammar.get(0);
        assert_eq!(result, &text);
    }
}
