use lazy_static::lazy_static;
use std::io::BufRead;

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
