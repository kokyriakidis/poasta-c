# poasta-c

C/C++ bindings for the [poasta](https://github.com/broadinstitute/poasta) library (Partial Order Alignment).

## Supported Features

- **Global Alignment**: Currently, only global alignment is supported. Semi-global (ends-free) alignment is not implemented in the upstream `poasta` library.
- **Two Gap Models**:
  - **Simple Affine Gap**: One gap open penalty and one gap extend penalty (`poasta_add_sequence` and `poasta_add_sequence_with_weight`)
  - **Two-Piece Affine Gap**: Two gap penalty pairs, automatically choosing the cheaper option for each gap (`poasta_add_sequence_2piece` and `poasta_add_sequence_2piece_with_weight`). Useful for better modeling of short vs long gaps.

## Building

This crate is designed to be built as a static or dynamic library.

```bash
cargo build --release
```

This will generate:
- `target/release/libpoasta_c.so` (Dynamic library)
- `target/release/libpoasta_c.a` (Static library)
- `poasta.h` (C header)

## Usage in C++

Include `poasta.h` and link against the library.

### Example

```cpp
#include "poasta.h"
#include <iostream>
#include <vector>
#include <cstring>

int main() {
    // 1. Create a new graph
    PoastaGraph* graph = poasta_create_graph();

    // 2. Add sequences
    
    // Simple Affine Gap Model
    // Parameters: graph, sequence, length, mismatch_score, gap_extend, gap_open
    // Default scoring: mismatch=4, gap_extend=2, gap_open=6
    // Default weight: 1 (each base gets weight 1)
    const char* seq1 = "ACGTACGT";
    poasta_add_sequence(graph, seq1, strlen(seq1), 4, 2, 6);

    const char* seq2 = "ACGTTCGT";
    poasta_add_sequence(graph, seq2, strlen(seq2), 4, 2, 6);

    // Simple affine with custom weight (useful for duplicate sequences)
    // If you have 5 identical sequences, you can add once with weight=5 instead of 5 calls
    // Default weight for non-weighted functions is 1
    const char* seq3 = "ACGTACGT";
    poasta_add_sequence_with_weight(graph, seq3, strlen(seq3), 5, 4, 2, 6);

    // Two-Piece Affine Gap Model (better for short vs long gaps)
    // Parameters: graph, sequence, length, mismatch_score, gap_extend1, gap_open1, gap_extend2, gap_open2
    // Default two-piece scoring: mismatch=4, gap_extend1=2, gap_open1=4, gap_extend2=1, gap_open2=24
    // Default weight: 1 (each base gets weight 1)
    const char* seq4 = "ACGTACGT";
    poasta_add_sequence_2piece(graph, seq4, strlen(seq4), 4, 2, 4, 1, 24);

    // Two-piece with custom weight
    // Default weight for non-weighted functions is 1
    const char* seq5 = "ACGTACGT";
    poasta_add_sequence_2piece_with_weight(graph, seq5, strlen(seq5), 3, 4, 2, 4, 1, 24);

    // 3. Get Multiple Sequence Alignment (MSA)
    PoastaMsa msa = poasta_get_msa(graph);
    
    std::cout << "MSA (" << msa.num_sequences << " sequences):" << std::endl;
    for (size_t i = 0; i < msa.num_sequences; ++i) {
        std::cout << msa.sequences[i] << std::endl;
    }
    
    // Free MSA memory
    poasta_free_msa(msa);

    // 4. Get GFA output (optional)
    char* gfa = poasta_get_gfa(graph);
    if (gfa) {
        std::cout << "\nGFA Output:\n" << gfa << std::endl;
        free(gfa); // Use standard free() for the string
    }

    // 5. Free the graph
    poasta_free_graph(graph);

    return 0;
}
```

## API Reference

### Types

- `PoastaGraph`: Opaque pointer to the POA graph.
- `PoastaMsa`: Struct containing the MSA result.
    - `sequences`: Array of C strings (`char**`).
    - `num_sequences`: Number of sequences.

### Functions

#### Graph Management

- `poasta_create_graph()`: Creates a new graph.
- `poasta_free_graph(graph)`: Frees the graph.

#### Simple Affine Gap Model

- `poasta_add_sequence(graph, seq, len, mismatch, gap_extend, gap_open)`: Adds a sequence using simple affine gap model (Global alignment). **Default weight: 1** (each base gets weight 1). This is equivalent to calling `poasta_add_sequence_with_weight` with `weight=1`.
- `poasta_add_sequence_with_weight(graph, seq, len, weight, mismatch, gap_extend, gap_open)`: Adds a sequence with a specified weight using simple affine gap model. This is useful when you have many identical sequences - instead of calling `poasta_add_sequence` multiple times, you can add the sequence once with a weight equal to the count of identical sequences. All bases in the sequence will have the same weight value.

#### Two-Piece Affine Gap Model

- `poasta_add_sequence_2piece(graph, seq, len, mismatch, gap_extend1, gap_open1, gap_extend2, gap_open2)`: Adds a sequence using two-piece affine gap model (Global alignment). This model uses two gap penalty pairs and automatically chooses the cheaper option for each gap. Better for modeling biological sequences where short and long gaps may have different characteristics. **Default weight: 1** (each base gets weight 1). This is equivalent to calling `poasta_add_sequence_2piece_with_weight` with `weight=1`.
- `poasta_add_sequence_2piece_with_weight(graph, seq, len, weight, mismatch, gap_extend1, gap_open1, gap_extend2, gap_open2)`: Adds a sequence with a specified weight using two-piece affine gap model. Same as above, but with a custom weight for the entire sequence.

**Two-Piece Gap Model Explanation**: The gap penalty for a gap of length ℓ is computed as `min(gap_open1 + ℓ × gap_extend1, gap_open2 + ℓ × gap_extend2)`. Typically, the first pair (gap_open1, gap_extend1) favors short gaps, while the second pair (gap_open2, gap_extend2) favors long gaps. For example, with `gap_open1=4, gap_extend1=2, gap_open2=24, gap_extend2=1`, short gaps use the first pair, while longer gaps switch to the second pair.

#### Output Functions

- `poasta_get_msa(graph)`: Generates the MSA. Caller must free result.
- `poasta_free_msa(msa)`: Frees the MSA memory.
- `poasta_get_gfa(graph)`: Returns GFA string. Caller must free result with `free()`.
