# poasta-c

C/C++ bindings for the [poasta](https://github.com/broadinstitute/poasta) library (Partial Order Alignment).

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
    // Parameters: graph, sequence, length, mismatch_score, gap_open, gap_extend
    // Default scoring: mismatch=4, gap_open=6, gap_extend=2
    
    const char* seq1 = "ACGTACGT";
    poasta_add_sequence(graph, seq1, strlen(seq1), 4, 6, 2);

    const char* seq2 = "ACGTTCGT";
    poasta_add_sequence(graph, seq2, strlen(seq2), 4, 6, 2);

    // 3. Add sequence with specific alignment mode (Global or SemiGlobal)
    const char* seq3 = "CGT";
    // SemiGlobal allows free gaps at ends of query and graph

    // 4. Get Multiple Sequence Alignment (MSA)
    PoastaMsa msa = poasta_get_msa(graph);
    
    std::cout << "MSA (" << msa.num_sequences << " sequences):" << std::endl;
    for (size_t i = 0; i < msa.num_sequences; ++i) {
        std::cout << msa.sequences[i] << std::endl;
    }
    
    // Free MSA memory
    poasta_free_msa(msa);

    // 5. Get GFA output (optional)
    char* gfa = poasta_get_gfa(graph);
    if (gfa) {
        std::cout << "\nGFA Output:\n" << gfa << std::endl;
        free(gfa); // Use standard free() for the string
    }

    // 6. Free the graph
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

- `poasta_create_graph()`: Creates a new graph.
- `poasta_free_graph(graph)`: Frees the graph.
- `poasta_add_sequence(graph, seq, len, mismatch, gap_open, gap_extend)`: Adds a sequence (Global alignment).
- `poasta_add_sequence_with_mode(...)`: Adds a sequence with a specific alignment mode.
- `poasta_get_msa(graph)`: Generates the MSA. Caller must free result.
- `poasta_free_msa(msa)`: Frees the MSA memory.
- `poasta_get_gfa(graph)`: Returns GFA string. Caller must free result with `free()`.
