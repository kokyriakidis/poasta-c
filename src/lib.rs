use std::ffi::CString;
use std::io::{Cursor, BufRead};
use std::os::raw::{c_char, c_int};
use std::slice;
use std::ptr;

use poasta::graphs::poa::POAGraph;
use poasta::aligner::scoring::{GapAffine, AlignmentType};
use poasta::aligner::config::AffineMinGapCost;
use poasta::aligner::PoastaAligner;
use poasta::io::fasta::poa_graph_to_fasta;
use poasta::io::graph::graph_to_gfa;

/// Opaque pointer to the POAGraph<u32>.
pub struct PoastaGraph(POAGraph<u32>);

/// Struct to hold the MSA result.
#[repr(C)]
pub struct PoastaMsa {
    pub sequences: *mut *mut c_char,
    pub num_sequences: usize,
}


/// Creates a new empty POAGraph.
#[unsafe(no_mangle)]
pub extern "C" fn poasta_create_graph() -> *mut PoastaGraph {
    let graph = POAGraph::<u32>::new();
    Box::into_raw(Box::new(PoastaGraph(graph)))
}

/// Frees the POAGraph.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_free_graph(graph: *mut PoastaGraph) {
    if !graph.is_null() {
        unsafe {
            drop(Box::from_raw(graph));
        }
    }
}

/// Adds a sequence to the graph (Global alignment).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_add_sequence(
    graph: *mut PoastaGraph,
    seq: *const c_char,
    len: usize,
    mismatch_score: u8,
    gap_open: u8,
    gap_extend: u8,
) -> c_int {
    if graph.is_null() || seq.is_null() {
        return -1;
    }

    let graph_inner = unsafe { &mut (*graph).0 };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    
    // Create a dummy name for the sequence (e.g. "seq_N")
    let seq_name = format!("seq_{}", graph_inner.sequences.len());
    let weights = vec![1; len];

    if graph_inner.is_empty() {
        // First sequence, just add it
        match graph_inner.add_alignment_with_weights(&seq_name, seq_slice, None, &weights) {
            Ok(_) => 0,
            Err(_) => -2,
        }
    } else {
        // Align and then add
        let scoring = GapAffine::new(mismatch_score, gap_extend, gap_open);
        
        // Always use Global alignment
        let aln_type = AlignmentType::Global;

        let aligner = PoastaAligner::new(AffineMinGapCost(scoring), aln_type);
        
        let result = aligner.align::<u32, _>(graph_inner, seq_slice);
        
        match graph_inner.add_alignment_with_weights(&seq_name, seq_slice, Some(&result.alignment), &weights) {
            Ok(_) => 0,
            Err(_) => -3,
        }
    }
}

/// Adds a sequence to the graph with a specified weight (Global alignment).
/// The weight applies to the entire sequence, meaning all bases will have the same weight.
/// This is useful when many identical sequences exist - instead of adding them multiple times,
/// you can add once with a weight representing the count.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_add_sequence_with_weight(
    graph: *mut PoastaGraph,
    seq: *const c_char,
    len: usize,
    weight: u32,
    mismatch_score: u8,
    gap_open: u8,
    gap_extend: u8,
) -> c_int {
    if graph.is_null() || seq.is_null() {
        return -1;
    }

    let graph_inner = unsafe { &mut (*graph).0 };
    let seq_slice = unsafe { slice::from_raw_parts(seq as *const u8, len) };
    
    // Create a dummy name for the sequence (e.g. "seq_N")
    let seq_name = format!("seq_{}", graph_inner.sequences.len());
    // Use the provided weight for all bases in the sequence
    let weights = vec![weight as usize; len];

    if graph_inner.is_empty() {
        // First sequence, just add it
        match graph_inner.add_alignment_with_weights(&seq_name, seq_slice, None, &weights) {
            Ok(_) => 0,
            Err(_) => -2,
        }
    } else {
        // Align and then add
        let scoring = GapAffine::new(mismatch_score, gap_extend, gap_open);
        
        // Always use Global alignment
        let aln_type = AlignmentType::Global;

        let aligner = PoastaAligner::new(AffineMinGapCost(scoring), aln_type);
        
        let result = aligner.align::<u32, _>(graph_inner, seq_slice);
        
        match graph_inner.add_alignment_with_weights(&seq_name, seq_slice, Some(&result.alignment), &weights) {
            Ok(_) => 0,
            Err(_) => -3,
        }
    }
}

/// Generates the MSA from the graph.
/// Returns a PoastaMsa struct. Caller must free it with poasta_free_msa.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_get_msa(graph: *mut PoastaGraph) -> PoastaMsa {
    if graph.is_null() {
        return PoastaMsa { sequences: ptr::null_mut(), num_sequences: 0 };
    }

    let graph_inner = unsafe { &(*graph).0 };
    let mut buffer = Vec::new();
    
    // Write MSA to buffer in FASTA format
    if let Err(_) = poa_graph_to_fasta(graph_inner, &mut buffer) {
        return PoastaMsa { sequences: ptr::null_mut(), num_sequences: 0 };
    }

    // Parse FASTA buffer to extract sequences
    let cursor = Cursor::new(buffer);
    let mut reader = std::io::BufReader::new(cursor);
    let mut line = String::new();
    let mut sequences = Vec::new();
    let mut current_seq = String::new();
    let mut in_seq = false;

    // Simple FASTA parser
    while reader.read_line(&mut line).unwrap() > 0 {
        let trimmed = line.trim();
        if trimmed.starts_with('>') {
            if in_seq {
                sequences.push(current_seq.clone());
                current_seq.clear();
            }
            in_seq = true;
        } else {
            if in_seq {
                current_seq.push_str(trimmed);
            }
        }
        line.clear();
    }
    if in_seq && !current_seq.is_empty() {
        sequences.push(current_seq);
    }

    // Convert to C strings
    let mut c_seqs = Vec::with_capacity(sequences.len());
    for s in sequences {
        let c_str = CString::new(s).unwrap();
        c_seqs.push(c_str.into_raw());
    }

    let ptr = c_seqs.as_mut_ptr();
    let len = c_seqs.len();
    std::mem::forget(c_seqs);

    PoastaMsa {
        sequences: ptr,
        num_sequences: len,
    }
}

/// Returns the graph in GFA format as a C string.
/// The caller must free the string using free().
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_get_gfa(graph: *mut PoastaGraph) -> *mut c_char {
    if graph.is_null() {
        return ptr::null_mut();
    }

    let graph_inner = unsafe { &(*graph).0 };
    let mut buffer = Vec::new();

    if let Err(_) = graph_to_gfa(&mut buffer, graph_inner) {
        return ptr::null_mut();
    }

    let s = String::from_utf8(buffer).unwrap_or_default();
    let c_str = CString::new(s).unwrap();
    c_str.into_raw()
}

/// Frees the PoastaMsa.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poasta_free_msa(msa: PoastaMsa) {
    if !msa.sequences.is_null() {
        let slice = unsafe { slice::from_raw_parts_mut(msa.sequences, msa.num_sequences) };
        for &mut ptr in slice {
            if !ptr.is_null() {
                unsafe {
                    drop(CString::from_raw(ptr));
                }
            }
        }
        unsafe {
            drop(Vec::from_raw_parts(msa.sequences, msa.num_sequences, msa.num_sequences));
        }
    }
}
