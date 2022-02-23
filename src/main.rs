use seq_io::fasta::{Reader, write_to};
use std::env::args;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use indexmap::IndexMap;

const EC_K: usize = 10;
const EC_CT: usize = 6;
const DIRS: [Dir; 2] = [Dir::Fwd, Dir::Bwd];
const GRAPH_K: usize = 40;
const CONTIG_TRESH: usize = 300;
const TIP_TRESH: usize = 80;

fn read_fasta(path: &str) -> Vec<Vec<u8>>{
    let mut reader = Reader::from_path(path).unwrap();
    let mut res = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        res.push(record.full_seq().into_owned());
    }
    res
}

struct Histogram {
    data: HashMap<Vec<u8>, usize>,
}
impl Histogram {
    fn new (reads: &Vec<Vec<u8>>) -> Self {
        let mut data = HashMap::new();
        for read in reads {
            for kmer in read.windows(EC_K) {
                let key = kmer.to_vec();
                let counter = data.entry(key).or_insert(0);
                *counter += 1;
            }
        }
        Histogram {data: data} 
    }
    fn get (&self, kmer: &[u8]) -> usize {
        self.data.get(kmer).cloned().unwrap_or(0)
    }
}

struct Graph {
    vertices: IndexMap<Vec<u8>, Vertex>,
}

#[derive(Default)]
struct Vertex {
    bwd: Vec<Edge>,
    fwd: Vec<Edge>,
}

#[derive(Copy, Clone)]
struct Edge {
    target_vertex: usize,
    weight: usize,
}

#[derive(Copy, Clone, PartialEq, Eq)]
enum Dir {
    Fwd,
    Bwd,
}

impl Dir {
    fn flip(self) -> Dir {
        match self {
            Dir::Fwd => Dir::Bwd,
            Dir::Bwd => Dir::Fwd,
        }
    }
}
impl Vertex {
    fn edges(&self, dir: Dir) -> &Vec<Edge> {
        match dir {
            Dir::Fwd => &self.fwd,
            Dir::Bwd => &self.bwd,
        }
    }
    fn edges_mut(&mut self, dir: Dir) -> &mut Vec<Edge> {
        match dir {
            Dir::Fwd => &mut self.fwd,
            Dir::Bwd => &mut self.bwd,
        }
    }
    fn add_edge(&mut self, dir: Dir, target_vertex: usize) {
        let edges = self.edges_mut(dir);
        for edge in edges.iter_mut() {
            if edge.target_vertex == target_vertex {
                edge.weight += 1;
                return;
            }
        }
        edges.push(Edge {
            target_vertex,
            weight: 1,
        });
    }
    fn remove_edge(&mut self, dir: Dir, target_vertex: usize) {
        self.edges_mut(dir).retain(|edge| edge.target_vertex != target_vertex);
    }

}
impl Graph {
    fn get_vertex(&mut self, kmer: &[u8]) -> usize {
        let e = self.vertices.entry(kmer.to_vec());
        let res = e.index();
        e.or_default();
        res
    }
    fn add_edge(&mut self, v1: usize, v2: usize) {
        self.vertices[v1].add_edge(Dir::Fwd, v2);
        self.vertices[v2].add_edge(Dir::Bwd, v1);
    }
    fn remove_vertex(&mut self, v: usize) {
        let vertex = std::mem::take(&mut self.vertices[v]);
        for dir in DIRS {
            for edge in vertex.edges(dir) {
                self.vertices[edge.target_vertex].remove_edge(dir.flip(), v);
            }
        }
    }
    fn get_string(&self, v: usize) -> &[u8] {
        self.vertices.get_index(v).unwrap().0
    }

    fn new(reads: &Vec<Vec<u8>>) -> Self {
        let mut res = Graph {vertices: IndexMap::new()};
        for read in reads {
            for pos in 0..read.len()-GRAPH_K {
                let v1 = res.get_vertex(&read[pos..pos+GRAPH_K]);
                let v2 = res.get_vertex(&read[pos+1..pos+1+GRAPH_K]);
                res.add_edge(v1, v2);
            }
        }
        res
    }
    fn assemble(&self) -> Vec<Vec<u8>>{
        let mut queue = Vec::new();
        for v in 0..self.vertices.len() {
            let vertex = &self.vertices[v];
            let weights_fwd: usize = vertex.fwd.iter().map(|e| e.weight).sum();
            let weights_bwd: usize = vertex.bwd.iter().map(|e| e.weight).sum();
            queue.push((weights_fwd + weights_bwd, v));
        }
        queue.sort_by(|a, b| b.0.cmp(&a.0));
        let mut processed = HashSet::new();
        let mut res = Vec::new();
        for (_, v) in queue {
            if processed.contains(&v) {
                continue;
            }
            processed.insert(v);
            let mut path = vec![v];
            for dir in [Dir::Bwd, Dir::Fwd] {
                let mut cv = v;
                loop {
                    let vertex = &self.vertices[cv];
                    match vertex.edges(dir).iter().filter(|edge| ! processed.contains(&edge.target_vertex)).max_by_key(|edge| edge.weight) {
                        None => break,
                        Some(edge) => {
                            path.push(edge.target_vertex);
                            processed.insert(edge.target_vertex);
                            cv = edge.target_vertex;
                        }
                    }
                }
                if dir == Dir::Bwd {
                    path.reverse();
                }
            }
            let mut contig = self.get_string(path[0]).to_vec();
            for &vertex in &path[1..] {
                let s = self.get_string(vertex);
                contig.push(*s.last().unwrap());
            }
            if contig.len() >= CONTIG_TRESH {
                res.push(contig);
            }
        }
        res
    }

    fn clean_tips(&mut self) {
        for v in 0..self.vertices.len() {
            for dir in DIRS {
                let vertex = &self.vertices[v];
                if ! vertex.edges(dir.flip()).is_empty() {
                    continue;
                }
                let mut path = vec![v];
                let mut cv = v;
                loop {
                    let curr_v = &self.vertices[cv];
                    let edges = curr_v.edges(dir);
                    if edges.len() != 1 {
                        break;
                    }
                    let tv = &self.vertices[edges[0].target_vertex];
                    let back_edges = tv.edges(dir.flip());
                    if back_edges.len() != 1 {
                        let best = back_edges.iter().map(|edge| edge.weight).max().unwrap();
                        // heuristic
                        if (edges[0].weight < best || edges[0].weight == 1) && path.len() < TIP_TRESH {
                            for v in path {
                                self.remove_vertex(v);
                            }
                        }
                        break;
                    }
                    path.push(edges[0].target_vertex);
                    cv = edges[0].target_vertex;
                }

            }
        }
    }

}

fn score_position(read: &[u8], pos: usize, hist: &Histogram) -> usize {
    let mut res = 1;
    for i in 0..EC_K {
        if i > pos {
            continue;
        }
        let pos_start = pos - i;
        let pos_end = pos_start + EC_K;
        if pos_end > read.len() {
            continue;
        }
        let kmer = &read[pos_start..pos_end];
        let freq = hist.get(kmer);
        res += freq.min(EC_CT);
    }
    res
}

fn correct_errors(reads: &mut Vec<Vec<u8>>){
    let kmer_hist = Histogram::new(reads);
    for read in reads {
        loop {
            let mut candidates = Vec::new();
            for pos in 0..read.len() {
                let orig_score = score_position(read, pos, &kmer_hist);
                let orig_val = read[pos];
                for &new_val in b"ACTGU" {
                    if new_val == orig_val {
                        continue;
                    }
                    read[pos] = new_val;
                    let new_score = score_position(read, pos, &kmer_hist);
                    // we want 50 percent better 
                    if 2 * new_score > 3 * orig_score {
                        candidates.push((new_score * 1000 / orig_score, pos, new_val));
                    }
                }
                read[pos] = orig_val;  
            }
            match candidates.into_iter().max() {
                Some((_, pos, new_val)) => {
                    read[pos] = new_val;
                }
                None => {
                    break;
                }
            }
        }
    }
}


fn main() {
    let path = args().nth(1).unwrap();
    let mut reads = read_fasta(&path);
    correct_errors(&mut reads);
    correct_errors(&mut reads);
    let mut graph = Graph::new(&reads);
    graph.clean_tips();
    let res = graph.assemble();
    let output_path = args().nth(2).unwrap();
    let mut output = File::create(output_path).unwrap();
    for (i, contig) in res.iter().enumerate() {
        write_to(&mut output, format!("CONTIG{}", i).as_bytes(), contig).unwrap();
    }
}
