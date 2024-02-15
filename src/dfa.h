#ifndef DFA_H_
#define DFA_H_

#include <vector>
#include <unordered_map>
#include <ctime>

#include "common.h"
#include "matrix.h"
#include "hashset.h"
#include "adjlist.h"
#include "edgeheap.h"
#include "edgearray.h"
#include "random.h"

const u32 kNoState = (u32)-1;

enum HashFunction {
    // Pick smallest state in a random permutation.
    MinHashOne,

    // Pick two smallest states in a random permutation, combine with bitwise-OR and modulo.
    MinHashTwo,

    // Pick k smallest states in a random permutation, combine with string hashing modulo.
    MinHashK
};

struct DFA {
    i32 alphabet_size;
    i32 state_count;

    // Transition table. Laid out sequentially, row by row.
    Matrix<state_t> txs;

    // A counter per state to indicate how many rules (if any) it accepts.
    i32 *accepting_states;
};

// NOTE: To be even more efficient we should layout these in a SOA fashion.
//       To determine the correct transition we just need to inspect the labels,
//       of which we can fit a lot in a cache line (64 bytes / 1 byte per char = 64 in line).
struct LabeledTx {
    state_t dest;
    label_t label;
};

struct Range {
    i32 index;
    i32 count;
};

struct D2FA {
    i32 alphabet_size;
    i32 state_count;

    // Labeled transitions. Stored as a big array where each state points to a contigious range.
    LabeledTx *txs = NULL;  // Full collection of all transitions.
    Range *tx_ranges = NULL; // txs[state] = range into txs that contain all the txs for state.

    // Default transitions. Either a state or kNoState for no default transitions.
    state_t *default_txs = NULL;

    // A counter per state to indicate how many rules (if any) it accepts.
    i32 *accepting_states = NULL;

    bool failed_to_construct = false; // Flag to indicate time-out.
    i32 stat_heap_resizes = 0; // Counter for how many times the edge heap resized during construction.
};

enum MinHashBoostApproach {
    KElements,         // Pick k smallest elements in one permutation.
    KPermutations,     // Pick smallest element in k permutations.
    KSamplesReplace,   // Sample k random edges with replacement.
    KSamplesNoReplace, // Sample k random edges without replacement.
};

struct D2FAParams {
    u32 path_length_bound = 0; // Zero = unbounded.

    bool do_use_sparse_srg = false; // False = use full SRG. In that case, below parameters are ignored.

    u32 min_hash_r = 0; // Number of rounds in SRG generation.
    u32 min_hash_k = 0; // Number of labels considered when hashing (per state).
    MinHashBoostApproach hash_approach = MinHashBoostApproach::KSamplesNoReplace; // Probability boosting approach.

};

void generate_primes();

void load_dfa(DFA *dfa, const char *path);
void save_dfa(const DFA &dfa, const char *path);

//void dfa_to_d2fa_with_rounds(const DFA &dfa, D2FA &out, const D2FAParams &params);
void dfa_to_d2fa_kruskal(const DFA &dfa, D2FA &out, const D2FAParams &params);
void dfa_to_d2fa_cutting(const DFA &dfa, D2FA &out, const D2FAParams &params);

void dfa_to_adfa_full_srg(const DFA &dfa, D2FA &out, const D2FAParams &params);
void dfa_to_adfa_sparse(const DFA &dfa, D2FA &out, const D2FAParams &params);

DFA merge(const DFA &a, const DFA &b);

void measure_bucket_sizes(DFA &dfa, HashFunction hash_type, u32 seed = 0, u32 min_hash_k = 1);

state_t match_character(const DFA &f,  state_t cur_state, label_t c);
state_t match_character(const D2FA &f, state_t cur_state, label_t c);

bool match(const DFA &dfa, const char *str);
bool match(const D2FA &f, const char *str);

void free_dfa(DFA &dfa);
void free_d2fa(D2FA &d2fa);

i32 count_common_txs(const DFA &dfa, i32 i, i32 j);

void hash_states_k_smallest(const DFA &dfa, u32 min_hash_k, i32 *out_hashes);
void hash_states_k_combined(const DFA &dfa, u32 min_hash_k, i32 *out_hashes);
void hash_states_k_samples_with_replace(const DFA &dfa, u32 min_hash_k, i32 *out_hashes);
void hash_states_k_samples_no_replace(const DFA &dfa, u32 min_hash_k, i32 *out_hashes);

// NOTE: Due to templates we need to place this function implementation here.
// This constructs the partial SRG.

// NOTE: Important that the endpoints of this edge are in different trees!
//       Otherwise the diameter calculation could overshoot.
const i32 kBigHeapWeight = 1 << 20;

struct SRGTimingInfo {
    clock_t clocks_to_hash;
    clock_t clocks_to_add_edges;
};

inline void clear_edge_container(EdgeHeap &out) { eheap_clear(out); }
inline void clear_edge_container(EdgeArray &out) { earray_clear(out); }
inline void clear_edge_container(AdjacencyList &out) { reset_adj_list(out); }
inline void clear_edge_container(Vector<Edge> &out) { out.clear(); }

inline void finalize_edge_container(AdjacencyList &out) { finalize_adj_list(out); }
inline void finalize_edge_container(Vector<Edge> &out) { /* Nothing to do */ }
inline void finalize_edge_container(EdgeHeap &out) { /* Nothing to do */ }
inline void finalize_edge_container(EdgeArray &out) { /* Nothing to do */ }

inline void add_edge(const DFA &dfa, AdjacencyList &out, i32 a, i32 b) {
    count_edge(out, a, b);
}
inline void add_edge(const DFA &dfa, Vector<Edge> &out, i32 a, i32 b) {
    i32 common_txs = count_common_txs(dfa, a, b);
    if (common_txs < 128) return;

    Edge edge;
    edge.source = a;
    edge.dest   = b;
    edge.heap_weight = common_txs;

    out.push_back(edge);
}
inline void add_edge(const DFA &dfa, EdgeHeap &out, i32 a, i32 b) {
    i32 common_txs = count_common_txs(dfa, a, b);
    if (common_txs < 128) return;

    SmallEdge edge;
    edge.source = a;
    edge.dest   = b;
    // edge.heap_weight = common_txs;

    i32 heap_weight = (1 + kBigHeapWeight) * common_txs; // NOTE: Corresponds to calculate_heap_weight for radius=0 (initial).
    eheap_push(out, edge, heap_weight);
}
inline void add_edge(const DFA &dfa, EdgeArray &out, i32 a, i32 b) {
    i32 common_txs = count_common_txs(dfa, a, b);
    if (common_txs < 128) return;

    Edge edge;
    edge.source = a;
    edge.dest   = b;
    edge.heap_weight = common_txs;

    earray_push(out, edge);
}

template <typename ContainerT>
i64 construct_sparse_srg_in_rounds(const DFA &dfa, u32 min_hash_k, MinHashBoostApproach boost_approach, u32 min_hash_rounds, ContainerT &out, SRGTimingInfo &timing_info) {
    i32 *hashes = alloc_array<i32>(dfa.state_count); // Hash for each state.

    Vector<i32> bucket_hashes; // All hash keys of non-empty buckets.
    UnorderedMap<i32, Vector<state_t>> buckets; // Map from hash to bucket contents.
    HashSet edge_set = hset_init(dfa.state_count); // Tracks which edges have already been added.

    i64 total_edges_for_srg = 0;

    clear_edge_container(out);

    timing_info = {}; // Zero.

    // printf("   Connecting initial graph.\n");

    // Ensure graph is connected.
    for (i32 i = 1; i < dfa.state_count; i++) {
        add_edge(dfa, out, 0, i);
        hset_insert(edge_set, 0, i);
        total_edges_for_srg++;
    }

    for (u32 r = 0; r < min_hash_rounds; r++) {
        // printf("   Start round %d / %d (edges=%ld, mbs=%ld):\n", r, min_hash_rounds - 1, total_edges_for_srg, mem_count.current / (1024*1024));

        clock_t t0 = clock();

        // printf("    Hashing states.\n");

        if      (boost_approach == MinHashBoostApproach::KElements)         hash_states_k_smallest(dfa, min_hash_k, hashes);
        else if (boost_approach == MinHashBoostApproach::KPermutations)     hash_states_k_combined(dfa, min_hash_k, hashes);
        else if (boost_approach == MinHashBoostApproach::KSamplesReplace)   hash_states_k_samples_with_replace(dfa, min_hash_k, hashes);
        else if (boost_approach == MinHashBoostApproach::KSamplesNoReplace) hash_states_k_samples_no_replace(dfa, min_hash_k, hashes);
        else assert(false);

        clock_t t1 = clock();

        // printf("    Constructing buckets.\n");

        // Construct the buckets.
        bucket_hashes.clear();
        buckets.clear();
        for (i32 i = 1; i < dfa.state_count; i++) {
            const i32 h = hashes[i];

            Vector<state_t> &bucket = buckets[h];
            if (bucket.empty()) bucket_hashes.push_back(h);
            bucket.push_back(i);
        }

        // printf("    Constructing edges.\n");

        // For each bucket, construct the edges.
        for (i32 b : bucket_hashes) {
            const Vector<state_t> &bucket = buckets[b];
            const u32 bucket_size = bucket.size();

            assert(bucket_size >= 1);
            if (bucket_size == 1) continue;

            // Sample a single edge per node.
            for (u32 i = 0; i < bucket_size; i++) {
                u32 j = (i + 1 + (random_u32() % (bucket_size - 1))) % bucket_size;
                if (hset_contains(edge_set, bucket[i], bucket[j])) continue;

                hset_insert(edge_set, bucket[i], bucket[j]);
                add_edge(dfa, out, bucket[i], bucket[j]);

                total_edges_for_srg++;
            }
        }

        clock_t t2 = clock();
        timing_info.clocks_to_hash      += t1 - t0;
        timing_info.clocks_to_add_edges += t2 - t1;

        f64 secs = (t2 - t1) / (f64)CLOCKS_PER_SEC;
        f64 ms = secs * 1000.0;
        f64 us = ms * 1000.0;
        f64 est_secs_left = (f64)((min_hash_rounds - r) * secs);
        // printf("    Elapsed: Total=%.2fms PerState=%.4fus (EstSecsLeft=%.2f)\n", ms, us / (f32)dfa.state_count, est_secs_left);
    }

    // printf("   Finalizing edge containers.\n");

    clock_t t1 = clock();
    finalize_edge_container(out);
    clock_t t2 = clock();
    timing_info.clocks_to_add_edges += t2 - t1;

    hset_free(edge_set);
    free_and_reset_array(hashes);

    assert(total_edges_for_srg == (i64)edge_set.count);
    return total_edges_for_srg;
}

template <typename ContainerT>
i64 construct_sparse_srg_in_rounds(const DFA &dfa, const D2FAParams &params, ContainerT &out, SRGTimingInfo &timing_info) {
    return construct_sparse_srg_in_rounds(dfa, params.min_hash_k, params.hash_approach, params.min_hash_r, out, timing_info);
}

#endif // DFA_H_
