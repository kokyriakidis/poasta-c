#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>

/// Opaque pointer to the POAGraph<u32>.
struct PoastaGraph;

/// Struct to hold the MSA result.
struct PoastaMsa {
  char **sequences;
  uintptr_t num_sequences;
};

extern "C" {

/// Creates a new empty POAGraph.
PoastaGraph *poasta_create_graph();

/// Frees the POAGraph.
void poasta_free_graph(PoastaGraph *graph);

/// Adds a sequence to the graph (Global alignment).
int poasta_add_sequence(PoastaGraph *graph,
                        const char *seq,
                        uintptr_t len,
                        uint8_t mismatch_score,
                        uint8_t gap_open,
                        uint8_t gap_extend);

/// Generates the MSA from the graph.
/// Returns a PoastaMsa struct. Caller must free it with poasta_free_msa.
PoastaMsa poasta_get_msa(PoastaGraph *graph);

/// Returns the graph in GFA format as a C string.
/// The caller must free the string using free().
char *poasta_get_gfa(PoastaGraph *graph);

/// Frees the PoastaMsa.
void poasta_free_msa(PoastaMsa msa);

}  // extern "C"
