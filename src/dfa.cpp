#include "dfa.h"

#include <string>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <ctime>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <deque>

#include "common.h"
#include "hash.h"
#include "matrix.h"
#include "util.h"
#include "unionfind.h"
#include "bitarray.h"
#include "adjlist.h"
#include "edgeheap.h"
#include "edgearray.h"
#include "hashset.h"
#include "log.h"
#include "random.h"

const i32 MIN_COMMON_TX_COUNT = 128;

struct PrimEdge {
    state_t source;
    state_t dest;

    i32 weight;
    i32 dist;
};

struct EdgeTrim {
    state_t source;
    state_t dest;
};

struct TreeBFSNode {
    state_t node;
    i32 parent;
};

// Node used in the traversal to construct defaults txs.
struct DefaultTxNode {
    state_t node;
};

struct OrderEdgeByWeightGreaterThan {
    inline bool operator() (const Edge &a, const Edge &b) {
        return a.heap_weight > b.heap_weight;

        // if (a.heap_weight != b.heap_weight) return a.heap_weight > b.heap_weight;
        // else if (a.source != b.source) return a.source < b.source;
        // else return a.dest < b.dest;
    }
};

i32 calculate_edge_weight(i32 a, i32 b, i32 common_txs, const i32 *radii) {
    i32 diameter = 1 + radii[a] + radii[b];
    return (1 + kBigHeapWeight) * common_txs - diameter;
}

inline i32 weigh_prim_edge(const PrimEdge &e) {
    return e.weight - (1 << e.dist);
}
struct OrderPrimHeuristic {
    inline bool operator() (const PrimEdge &a, const PrimEdge &b) {
        return weigh_prim_edge(a) < weigh_prim_edge(b);
    }
};
struct OrderPrimOnlyByWeight {
    inline bool operator() (const PrimEdge &a, const PrimEdge &b) {
        return a.weight < b.weight;
    }
};

bool eat_on_line(FILE *file, char target) {
    char c;
    do {
        c = fgetc(file);
        if (c == target) return true;
    } while (c == ' ');
    fprintf(stderr, "[ERROR] Did not find expected '%c' on line!\n", target);
    return false;
}

void load_dfa(DFA *dfa, const char *path) {
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "[ERROR] Failed to open DFA file '%s'\n", path);
        return;
    }

    // Skip comment line.
    char c;
    do {
        c = fgetc(file);
    } while (c != '\n');

    if (0 == fscanf(file, "%d", &dfa->state_count)) {
        fprintf(stderr, "[ERROR] Failed to read DFA state count in '%s'\n", path);
        return;
    }
    assert(dfa->state_count > 0);

    if (0 == fscanf(file, "%d", &dfa->alphabet_size)) {
        fprintf(stderr, "[ERROR] Failed to read DFA alphabet size in '%s'\n", path);
        return;
    }

    const i32 kMaxAlphabetSize = 1 << 8*sizeof(label_t);
    if (dfa->alphabet_size > kMaxAlphabetSize) {
        fprintf(stderr, "[ERROR] Alphabet size (%d) is larger than supported (%d)!\n", dfa->alphabet_size, kMaxAlphabetSize);
        return;
    }
    assert(dfa->alphabet_size > 0);

    // printf("  Loaded sizes for DFA: #states=%d alphabet=%d\n", dfa->state_count, dfa->alphabet_size);

    dfa->accepting_states = alloc_array<i32>(dfa->state_count);
    dfa->txs = alloc_matrix<state_t>(dfa->alphabet_size, dfa->state_count);

    for (i32 i = 0; i < dfa->state_count; i++) {
        fgetc(file); // Skip newline.

        if (!eat_on_line(file, '(')) {
            fprintf(stderr, "[ERROR] State %d: Expected outer '(' but didn't find it.\n", i);
            return;
        }
        if (!eat_on_line(file, '(')) {
            fprintf(stderr, "[ERROR] State %d: Expected inner '(' but didn't find it.\n", i);
            return;
        }

        i32 dest = 0;
        for (i32 j = 0; j < dfa->alphabet_size; j++) {
            if (0 == fscanf(file, "%d", &dest)) {
                fprintf(stderr, "[ERROR] State %d: Failed to read state for symbol %d (%c).\n", i, j, j);
                return;
            }
            at(dfa->txs, i, j) = dest;
        }
        if (!eat_on_line(file, ')')) {
            fprintf(stderr, "[ERROR] State %d: Expected inner ')' but didn't find it.\n", i);
            return;
        }

        i32 accepted_rule = 0;
        i32 accepting_rules_count = 0;
        while (0 < fscanf(file, "%d", &accepted_rule)) accepting_rules_count++;
        dfa->accepting_states[i] = accepting_rules_count;
        assert(accepting_rules_count >= 0);

        if (!eat_on_line(file, ')')) {
            fprintf(stderr, "[ERROR] State %d: Expected outer ')' but didn't find it.\n", i);
            return;
        }
    }

    assert_dfa_size(dfa);
    fclose(file);
}

void save_dfa(const DFA &dfa, const char *path) {
    FILE *file = fopen(path, "w");
    if (file == NULL) {
        fprintf(stderr, "[ERROR] Failed to create output DFA file '%s'\n", path);
        return;
    }

    fprintf(file, "# DFA dump\n");
    fprintf(file, "%d\n", dfa.state_count);
    fprintf(file, "%d\n", dfa.alphabet_size);

    for (i32 i = 0; i < dfa.state_count; i++) {
        fprintf(file, "( ");
        fprintf(file, " ( ");
        for (i32 j = 0; j < dfa.alphabet_size; j++) {
            fprintf(file, "%d ", at(dfa.txs, i, j));
        }
        fprintf(file, ")");

        // NOTE: We do not currently store the set of rules that an accepting state corresponds to.
        // Therefore we just store -1 for each state. We do this to match the DFA format.
        for (i32 j = 0; j < dfa.accepting_states[i]; j++) {
            fprintf(file, "-1 ");
        }
        fprintf(file, ")\n");
    }
    fclose(file);
}

void construct_default_txs_from_outside_going_in(D2FA &out, AdjacencyList &mst_edge_set, i32 mst_edge_count) {
    DefaultTxNode *forest_traversal_queue   = alloc_array<DefaultTxNode>(out.state_count); // Order to visit leaves in.
    i32 *forest_traversal_unmatched_degrees = alloc_array<i32>(out.state_count); // The number of unmatched neighbours in traversal algorithm.
    BitArray forest_traversal_matched       = ba_make(out.state_count);

    // We orient all default transitions towards the center of the spanning tree.
    // We work from the leaves and inwards, creating default txs as we go and keeping track
    // of which nodes become "leaves" as their neighbour get defaults txs.
    // The last node to be visited will be a center (could be multiple), and promoted to root of the tree.
    memzero(forest_traversal_unmatched_degrees, out.state_count);
    ba_zero(forest_traversal_matched);
    i32 next_in_queue = 0; // Pointer to the next element in the queue.
    for (i32 i = 0; i < out.state_count; i++) {
        i64 a, b;
        get_edge_range(mst_edge_set, i, a, b);
        const i64 range_size = b - a;

        forest_traversal_unmatched_degrees[i] = range_size;
        if (range_size == 1) {
            forest_traversal_queue[next_in_queue].node      = i;
            next_in_queue++;
        }
    }

    i32 default_txs_added = 0; // Used for verification. Delete for speed.
    for (i32 i = 0; i < next_in_queue; i++) {
        if (default_txs_added == mst_edge_count) break; // NOTE: This will always be the termination condition.

        const DefaultTxNode &cur = forest_traversal_queue[i];
        i64 a, b;
        get_edge_range(mst_edge_set, cur.node, a, b);

        // The root of a tree can end up being queued when second-to-last child becomes matched,
        // but then when the last child gets matched the root is no longer a leaf, since it
        // has zero (not one) unmatched neighbours. In this case we simply skip it, as that
        // spanning tree is fully constructed then.
        if (forest_traversal_unmatched_degrees[cur.node] == 0) continue;

        assert(!ba_get(forest_traversal_matched, cur.node));
        assert(forest_traversal_unmatched_degrees[cur.node] == 1); // Should be "leaf".
        ba_set(forest_traversal_matched, cur.node, true); // Maybe easier to understand if we set this when tx is created?

        // Search to neighbours.
        i32 unmatched_neighbour_count = 0;
        for (i64 j = a; j < b; j++) {
            const i32 neighbour = mst_edge_set.dests[j];
            const bool neighbour_matched = ba_get(forest_traversal_matched, neighbour);
            if (neighbour_matched) continue;
            unmatched_neighbour_count++;

            // Add default transition towards (only) unmatched neighbour.
            // NOTE: We could break out of the inner loop here, as only one tx should be added.
            assert(out.default_txs[cur.node] == kNoState); // Should be a new default tx we construct.
            out.default_txs[cur.node] = neighbour; // Source / dest are bucket local indices.
            default_txs_added++;

            // Update neighbour.
            forest_traversal_unmatched_degrees[neighbour]--;
            if (forest_traversal_unmatched_degrees[neighbour] == 1) {
                forest_traversal_queue[next_in_queue].node = neighbour;
                next_in_queue++;
            } else assert(forest_traversal_unmatched_degrees[neighbour] >= 0); // Remove for speed.
        }

        assert(unmatched_neighbour_count == 1); // Otherwise we didn't process a leaf.
    }

    // NOTE: If we are operating on a fully connected SRG, then this check should be true.
    // assert(next_in_queue == out.state_count); // We have visited each node once.
    assert(default_txs_added == mst_edge_count); // Every edge should become a default tx.

    ba_delete(forest_traversal_matched);
    free_and_reset_array(forest_traversal_unmatched_degrees);
    free_and_reset_array(forest_traversal_queue);
}

i32 count_common_txs(const DFA &dfa, i32 i, i32 j) {
    i32 common_txs = 0;
    for (i32 l = 0; l < dfa.alphabet_size; l++) {
        common_txs += at(dfa.txs, i, l) == at(dfa.txs, j, l);
    }
    return common_txs;
}

// Copy in any labelled txs that have not been compressed.
// Requires that you know how many txs have been compressed (compressed_tx_count).
void copy_in_uncompressed_labelled_txs(const DFA &dfa, D2FA &out, i32 compressed_tx_count) {
    assert(compressed_tx_count >= 0);

    // Copy in the transitions.
    const i64 total_dfa_txs     = out.state_count * out.alphabet_size;
    const i64 total_labeled_txs = total_dfa_txs - compressed_tx_count;
    assert(total_labeled_txs >= 0);

    out.txs = alloc_array<LabeledTx>(total_labeled_txs);

    i64 tx_ptr = 0;
    i64 range_ptr = 0;
    for (i32 i = 0; i < out.state_count; i++) {
        const u32 default_tx_dest = out.default_txs[i];

        Range &range = out.tx_ranges[range_ptr];
        range.index = tx_ptr;

        if (default_tx_dest == kNoState) {
            for (i32 j = 0; j < out.alphabet_size; j++) {
                out.txs[tx_ptr].label = j;
                out.txs[tx_ptr].dest  = at(dfa.txs, i, j);
                tx_ptr++;
            }
        } else {
            for (i32 j = 0; j < out.alphabet_size; j++) {
                const state_t dest = at(dfa.txs, i, j);
                if (dest == at(dfa.txs, default_tx_dest, j)) continue;
                out.txs[tx_ptr].label = j;
                out.txs[tx_ptr].dest  = dest;
                tx_ptr++;
            }
        }

        range.count = tx_ptr - range.index;
        range_ptr++;
    };

    if(total_labeled_txs != tx_ptr) {
        printf(" ERROR: Total=%ld TxPtr=%ld\n", total_labeled_txs, tx_ptr);
        assert(false);
    }
}

void copy_in_uncompressed_labelled_txs_with_counting(const DFA &dfa, D2FA &out) {
    i64 compressed_tx_count = 0;
    for (i32 i = 0; i < out.state_count; i++) {
        const u32 default_tx_dest = out.default_txs[i];
        if (default_tx_dest == kNoState) continue;

        compressed_tx_count += count_common_txs(dfa, i, default_tx_dest);
    }

    const i64 total_dfa_txs     = out.state_count * out.alphabet_size;
    const i64 total_labeled_txs = total_dfa_txs - compressed_tx_count;
    assert(total_labeled_txs >= 0);
    out.txs = alloc_array<LabeledTx>(total_labeled_txs);

    i64 tx_ptr = 0;
    i64 range_ptr = 0;
    for (i32 i = 0; i < out.state_count; i++) {
        const u32 default_tx_dest = out.default_txs[i];

        Range &range = out.tx_ranges[range_ptr];
        range.index = tx_ptr;

        if (default_tx_dest == kNoState) {
            for (i32 j = 0; j < out.alphabet_size; j++) {
                out.txs[tx_ptr].label = j;
                out.txs[tx_ptr].dest  = at(dfa.txs, i, j);
                tx_ptr++;
            }
        } else {
            for (i32 j = 0; j < out.alphabet_size; j++) {
                const state_t dest = at(dfa.txs, i, j);
                if (dest == at(dfa.txs, default_tx_dest, j)) continue; // This must be slow...
                out.txs[tx_ptr].label = j;
                out.txs[tx_ptr].dest  = dest;
                tx_ptr++;
            }
        }

        range.count = tx_ptr - range.index;
        range_ptr++;
    };

    if(total_labeled_txs != tx_ptr) {
        printf(" ERROR: Total=%ld TxPtr=%ld\n", total_labeled_txs, tx_ptr);
        assert(false);
    }
}

// Hash by combining the k-smallest elements a single permutation.
void hash_states_k_smallest(const DFA &dfa, u32 min_hash_k, i32 *out_hashes) {
    u32 *mins = alloc_array<u32>(min_hash_k); // Buffer used for finding the k smallest (after permutation) neighbours of a state.

    const u64 a = uniform_u64_for_hash();
    const u64 b = uniform_u64_for_hash();

    // NOTE: For min_hash_k == 1, this is identical to k_combined.
    if (min_hash_k == 1) {
        for (i32 i = 0; i < dfa.state_count; i++) {
            out_hashes[i] = calculate_single_min_hash(dfa, i, a, b) % dfa.state_count;
        }
    } else {
        VectorHash vhash = vhash_init(min_hash_k);
        for (i32 i = 0; i < dfa.state_count; i++) {
            vhash_clear(vhash);

            calculate_k_min_hashes(dfa, i, a, b, min_hash_k, mins);
            for (u32 k = 0; k < min_hash_k / 2; k++) {
                vhash_add_pair(vhash, mins[2*k], mins[2*k+1]);
            }
            if (min_hash_k % 2 == 1) {
                vhash_add_single(vhash, mins[min_hash_k - 1]);
            }
            out_hashes[i] = vhash_finalize(vhash);
        }
        vhash_free(vhash);
    }

    free_and_reset_array(mins);
}

// Hash by combining the smallest element in k permutations.
void hash_states_k_combined(const DFA &dfa, u32 min_hash_k, i32 *out_hashes) {
    u64 *as = alloc_array<u64>(min_hash_k);
    u64 *bs = alloc_array<u64>(min_hash_k);
    for (u32 i = 0; i < min_hash_k; i++) {
        as[i] = uniform_u64_for_hash();
        bs[i] = uniform_u64_for_hash();
    }

    VectorHash vhash = vhash_init(min_hash_k);
    for (i32 i = 0; i < dfa.state_count; i++) {
        vhash_clear(vhash);
        for (u32 k = 0; k < min_hash_k / 2; k++) {
            const u32 h1 = calculate_single_min_hash(dfa, i, as[2*k],   bs[2*k]);
            const u32 h2 = calculate_single_min_hash(dfa, i, as[2*k+1], bs[2*k+1]);
            vhash_add_pair(vhash, h1, h2);
        }
        if (min_hash_k % 2 == 1) {
            const u32 h = calculate_single_min_hash(dfa, i, as[min_hash_k - 1], bs[min_hash_k - 1]);
            vhash_add_single(vhash, h);
        }
        out_hashes[i] = vhash_finalize(vhash);
    }
    vhash_free(vhash);

    free_and_reset_array(as);
    free_and_reset_array(bs);
}

// Hash by sampling k (not guaranteed to be unique) indices.
void hash_states_k_samples_with_replace(const DFA &dfa, u32 min_hash_k, i32 *out_hashes) {
    // These are the edge labels we sample.
    u32 *cs = alloc_array<u32>(min_hash_k);
    for (u32 i = 0; i < min_hash_k; i++) cs[i] = random_u32() % dfa.alphabet_size;

    if (min_hash_k == 1) {
        for (i32 i = 0; i < dfa.state_count; i++) {
            out_hashes[i] = at(dfa.txs, i , cs[0]);
        }
    } else {
        VectorHash vhash = vhash_init(min_hash_k);
        for (i32 i = 0; i < dfa.state_count; i++) {
            vhash_clear(vhash);
            for (u32 k = 0; k < min_hash_k / 2; k++) {
                vhash_add_pair(vhash, at(dfa.txs, i, cs[2*k]), at(dfa.txs, i, cs[2*k+1]));
            }
            if (min_hash_k % 2 == 1) {
                i32 x = at(dfa.txs, i, cs[min_hash_k - 1]);
                vhash_add_single(vhash, x);
            }
            out_hashes[i] = vhash_finalize(vhash);
        }
        vhash_free(vhash);
    }


    free_and_reset_array(cs);
}

// Hash by sampling k unique indices.
void hash_states_k_samples_no_replace(const DFA &dfa, u32 min_hash_k, i32 *out_hashes) {
    assert((i32)min_hash_k <= dfa.alphabet_size); // Without replacement we run out.

    // Tracks which symbols we already sampled.
    BitArray sampled = ba_make(dfa.alphabet_size);
    ba_zero(sampled);

    // These are the edge labels we sample.
    u32 *cs = alloc_array<u32>(min_hash_k);
    for (u32 i = 0; i < min_hash_k; i++) {
        u32 c = 0;
        do {
            c = random_u32() % dfa.alphabet_size;
        } while (ba_get(sampled, c));
        ba_set(sampled, c, true);
        cs[i] = c;
    }
    ba_delete(sampled);

    VectorHash vhash = vhash_init(min_hash_k);
    for (i32 i = 0; i < dfa.state_count; i++) {
        vhash_clear(vhash);
        for (u32 k = 0; k < min_hash_k / 2; k++) {
            // NOTE: If we use the same labels for every state, and a vector hashing algorithm (order matters),
            //       then there is no need to pack the label into the hash as it is implied by the order.
            vhash_add_pair(vhash, at(dfa.txs, i, cs[2*k]), at(dfa.txs, i, cs[2*k+1]));
        }
        if (min_hash_k % 2 == 1) {
            i32 x = at(dfa.txs, i, cs[min_hash_k - 1]);
            vhash_add_single(vhash, x);
        }
        out_hashes[i] = vhash_finalize(vhash);
    }
    vhash_free(vhash);

    free_and_reset_array(cs);
}

void update_radii(i32 source, i32 dest, i32 source_radius, i32 tree_size, i32 *radii, Vector<Vector<i32>> &forest_edges, EdgeHeap &edge_heap, const DFA &dfa, Vector<Vector<std::pair<i32,i32>>> &common_txs_cache) {
    i32 queue_start = 0;
    Vector<std::tuple<i32, i32>> queue; // (Node, Dist)
    UnorderedSet<i32> seen;

    queue.reserve(tree_size);
    seen.reserve(tree_size);

    // NOTE: radii[source] could have just changed, in a prior call, which is why we pass it.
    queue.push_back({ dest, 1 + source_radius });
    seen.insert(source); // We do not want to go "back" over the edge.
    seen.insert(dest);

    while (queue_start < (i32)queue.size()) {
        auto front = queue[queue_start];
        i32 cur                 = std::get<0>(front);
        i32 radius_through_root = std::get<1>(front);
        queue_start++;

        // If the radius has been reduced we update all incident edges in the heap.
        if (radii[cur] < radius_through_root) {
            radii[cur] = radius_through_root;

            if (common_txs_cache[cur].empty()) {
                for (i32 i = 0; i < dfa.state_count; i++) {
                    if (i == cur) continue;
                    u64 edge_key = to_edge_key(cur, i);
                    if (!eheap_contains(edge_heap, edge_key)) continue;
                    common_txs_cache[cur].push_back({ i, count_common_txs(dfa, cur, i) });
                }
            }

            for (auto next : common_txs_cache[cur]) {
                const i32 i          = next.first;
                const i32 common_txs = next.second;
                u64 edge_key = to_edge_key(cur, i);
                i32 new_weight = calculate_edge_weight(cur, i, common_txs, radii);
                eheap_update_weight(edge_heap, edge_key, new_weight);
            }
        }

        for (auto next : forest_edges[cur]) {
            if (seen.count(next)) continue;
            seen.insert(next);
            queue.push_back({ next, radius_through_root + 1 });
        }
    }
}

void dfa_to_d2fa_kruskal(const DFA &dfa, D2FA &out, const D2FAParams &params) {
    const i32 kMaxAlphabetSize = 1 << 8*sizeof(label_t);
    assert(dfa.alphabet_size <= kMaxAlphabetSize);

    out.alphabet_size = dfa.alphabet_size;
    out.state_count   = dfa.state_count;

    // NOTE: At this point we do not know the number of labelled transitions, so we do not allocate yet.
    // NOTE: We perform copies to construct a new D2FA without modifying the DFA.
    //       It might be more efficient to perform a transformation.
    out.tx_ranges        = alloc_array<Range>(out.state_count);
    out.default_txs      = alloc_array<state_t>(out.state_count);
    out.accepting_states = alloc_array<i32>(out.state_count);

    reset_mem_count();

    fill_array(out.default_txs, out.state_count, kNoState);
    copy_array(out.accepting_states, dfa.accepting_states, out.state_count);

    const clock_t t0 = clock();

    // We counting sort the edges as the weight universe is small.
    const i32 max_different_distances_count = out.alphabet_size + 1; // Number of different distance values.

    // We use UnionFind for Kruskals.
    UnionFind uf = uf_init(out.state_count);

    AdjacencyList tree_edge_set = make_adj_list(out.state_count);

    const clock_t t1 = clock();

    i64 tree_edge_count = 0;

    if (params.path_length_bound > 0) {
        i32 *radii = alloc_array<i32>(out.state_count); // Distance from the node to the furthest node in the same tree.

        SRGTimingInfo srg_time_info;
        EdgeHeap edge_heap = eheap_init(out.state_count);

        Vector<Vector<i32>> forest_edges(out.state_count);
        Vector<Vector<std::pair<i32,i32>>> common_txs_cache(out.state_count);

        // Construct the edges.
        memzero(radii, out.state_count);

        i64 total_edges_in_srg = 0;
        if (params.do_use_sparse_srg) {
            // printf("  Constructing sparse SRG...\n");
            construct_sparse_srg_in_rounds(dfa, params, edge_heap, srg_time_info);
            total_edges_in_srg = edge_heap.edges.size();
        } else {
            // printf("  Constructing full SRG...\n");
            for (i32 i = 0; i < out.state_count - 1; i++) {
                for (i32 j = i + 1; j < out.state_count; j++) {
                    // NOTE: This distance is not the amount actually saved, as we add a default tx.
                    // NOTE: Re-counting here is rather wasteful.
                    i32 common_txs = count_common_txs(dfa, i, j);
                    if (common_txs < MIN_COMMON_TX_COUNT) continue;

                    i32 distance = out.alphabet_size - common_txs;
                    assert(0 <= distance && distance < max_different_distances_count); // Remove for speed.

                    SmallEdge edge;
                    edge.source = i; // NOTE: Edge endpoints are indices into the bucket.
                    edge.dest   = j;
                    //edge.heap_weight = common_txs;

                    i32 heap_weight = calculate_edge_weight(i, j, common_txs, radii);
                    eheap_push(edge_heap, edge, heap_weight);

                    total_edges_in_srg++;
                }
            }
            //assert(srg_edge_count == total_edges_in_srg);
        }
        out.stat_heap_resizes = edge_heap.resize_count;

        // Run heuristic Kruskals.
        const i32 max_diameter = 2 * params.path_length_bound;

        // printf("  Running Kruskals...\n");
        while (!eheap_is_empty(edge_heap)) {
            // Pop max-weight edge.
            SmallEdge edge = eheap_pop(edge_heap);

            // Check if edge would cause a cycle.
            // We perform Find by iterating to the root of the tree.
            i32 id1 = uf_find(uf, edge.source);
            i32 id2 = uf_find(uf, edge.dest);
            if (id1 == id2) continue; // Same component, so it would form a cycle.

            // Check if we go above the diameter bound.
            i32 diameter = 1 + radii[edge.source] + radii[edge.dest];
            if (diameter > max_diameter) {
                //printf(" >>> Skip edge (%d - %d): Diam=1+%d+%d = %d\n", edge.source, edge.dest, radii[edge.source], radii[edge.dest], diameter);
                continue;
            }

            // Link in tree.
            forest_edges[edge.source].push_back(edge.dest);
            forest_edges[edge.dest].push_back(edge.source);

            // Update radii and weights for elements.
            i32 radius_source = radii[edge.source];
            i32 radius_dest   = radii[edge.dest];
            update_radii(edge.source, edge.dest,   radius_source, uf.szs[id2], radii, forest_edges, edge_heap, dfa, common_txs_cache);
            update_radii(edge.dest,   edge.source, radius_dest,   uf.szs[id1], radii, forest_edges, edge_heap, dfa, common_txs_cache);

            // Add to spanning tree.
            count_edge(tree_edge_set, edge.source, edge.dest);
            tree_edge_count++;

            uf_union(uf, id1, id2);
        }

        free_and_reset_array(radii);
    } else { // === Use simple non-heap based approach ===
        Vector<Edge> edges;
        edges.reserve(out.state_count);

        SRGTimingInfo srg_time_info;

        if (params.do_use_sparse_srg) {
            // printf("  Constructing sparse SRG...\n");
            construct_sparse_srg_in_rounds(dfa, params, edges, srg_time_info);
        } else {
            // printf("  Constructing full SRG...\n");
            i64 total_edges_in_srg = 0;
            for (i32 i = 0; i < out.state_count - 1; i++) {
                for (i32 j = i + 1; j < out.state_count; j++) {
                    i32 common_txs = count_common_txs(dfa, i, j);
                    if (common_txs < MIN_COMMON_TX_COUNT) continue;

                    Edge edge;
                    edge.source      = i;
                    edge.dest        = j;
                    edge.heap_weight = common_txs;
                    edges.push_back(edge);

                    total_edges_in_srg++;
                }
            }
        }

        // printf("  Sorting edges...\n");
        std::sort(edges.begin(), edges.end(), OrderEdgeByWeightGreaterThan());

        // printf("  Running Kruskals...\n");
        for (const Edge &edge : edges) {
            // Check if edge would cause a cycle.
            // We perform Find by iterating to the root of the tree.
            i32 id1 = uf_find(uf, edge.source);
            i32 id2 = uf_find(uf, edge.dest);
            if (id1 == id2) continue; // Same component, so it would form a cycle.

            // Add to spanning tree.
            count_edge(tree_edge_set, edge.source, edge.dest);
            tree_edge_count++;

            if (tree_edge_count == out.state_count - 1) break; // Tree fully constructed.

            uf_union(uf, id1, id2);
        }
    }

    // printf("  Finalizing adjacency list...\n");
    finalize_adj_list(tree_edge_set);

    // printf("  Constructing default txs...\n");
    construct_default_txs_from_outside_going_in(out, tree_edge_set, tree_edge_count);

    const clock_t t2 = clock();

    // printf("  Copying unlabelled txs...\n");
    copy_in_uncompressed_labelled_txs_with_counting(dfa, out);

    const clock_t t3 = clock();

    const f64 t_secs_total     = (f64)(t3 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_srg    = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_forest = (f64)(t2 - t1) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_txs    = (f64)(t3 - t2) / (f64)CLOCKS_PER_SEC;

    log(Log::Stat, "Time (secs): Total=%.4f\n\tSRGs=%.4f (%.2f%%)\n\tForest=%.4f (%.2f%%)\n\tTxs =%.4f (%.2f%%)\n",
           t_secs_total,
           t_secs_to_srg,    100.0 * t_secs_to_srg    / t_secs_total,
           t_secs_to_forest, 100.0 * t_secs_to_forest / t_secs_total,
           t_secs_to_txs,    100.0 * t_secs_to_txs    / t_secs_total);

    uf_delete(uf);
    free_ajd_list(tree_edge_set);
}

template <typename QueueT>
void prim_add_edges_from_node(i32 node, i32 node_dist, const AdjacencyList &neighbours, const DFA &dfa, QueueT &queue, BitArray &nodes_in_tree) {
    i64 a, b;
    get_edge_range(neighbours, node, a, b);
    for (i64 i = a; i < b; i++) {
        u32 dest = neighbours.dests[i];
        if (ba_get(nodes_in_tree, dest)) continue;

        // NOTE: There is an additional unnecesarry copy here, because we do not use emplace.
        i32 common_txs = count_common_txs(dfa, node, dest);

        PrimEdge e;
        e.source = node;
        e.dest   = dest;
        e.weight = common_txs;
        e.dist   = node_dist + 1;
        queue.push(e);
    }
}

void construct_preorder(u32 root, i32 &cur_index, i32 *out_order, const u32 *default_txs, const AdjacencyList &tree_edge_set) {
    out_order[cur_index++] = root;

    i64 a, b;
    get_edge_range(tree_edge_set, root, a, b);
    for (i64 i = a; i < b; i++) {
        const u32 dest = tree_edge_set.dests[i];
        if (default_txs[root] == dest) continue; // Do not go up in the tree.
        construct_preorder(dest, cur_index, out_order, default_txs, tree_edge_set);
    }
}

i32 bfs_to_exhaustion(TreeBFSNode *queue, i32 state_count, i32 root, AdjacencyList &alist) {
    BitArray seen = ba_make(state_count);
    ba_zero(seen);
    ba_set(seen, root, true);

    i32 i0 = 0;
    i32 i1 = 1;
    queue[0].node   = root;
    queue[0].parent = -1; // Index of parent in queue.

    while (i0 < i1) {
        const TreeBFSNode &cur = queue[i0];
        i64 a, b;
        get_edge_range(alist, cur.node, a, b);
        for (i64 i = a; i < b; i++) {
            u32 dest = alist.dests[i];
            if (ba_get(seen, dest)) continue; // Don't go back in tree.
            ba_set(seen, dest, true);

            queue[i1].node   = dest;
            queue[i1].parent = i0;
            i1++;
        }
        i0++;
    }
    assert(i1 == i0 && i0 <= state_count); // i0 == state_count only if single connected tree.

    ba_delete(seen);
    return i0 - 1;
}

void dfa_to_d2fa_cutting(const DFA &dfa, D2FA &out, const D2FAParams &params) {
    const i32 kMaxAlphabetSize = 1 << 8*sizeof(label_t);
    assert(dfa.alphabet_size <= kMaxAlphabetSize);

    out.alphabet_size = dfa.alphabet_size;
    out.state_count   = dfa.state_count;

    // NOTE: At this point we do not know the number of labelled transitions, so we do not allocate yet.
    // NOTE: We perform copies to construct a new D2FA without modifying the DFA.
    //       It might be more efficient to perform a transformation.
    out.tx_ranges        = alloc_array<Range>(out.state_count);
    out.default_txs      = alloc_array<state_t>(out.state_count);
    out.accepting_states = alloc_array<i32>(out.state_count);

    reset_mem_count();

    for (i32 i = 0; i < out.state_count; i++) out.default_txs[i] = kNoState;
    memcpy(out.accepting_states, dfa.accepting_states, sizeof(i32) * out.state_count);

    const clock_t t0 = clock();

    const i32 diameter_bound = params.path_length_bound > 0 ? 2 * params.path_length_bound : 0;

    SRGTimingInfo srg_time_info;

    // Construct the edges of the SRG.
    AdjacencyList srg_edge_set = make_adj_list(out.state_count);
    if (!params.do_use_sparse_srg) {
        // printf("  Constructing full SRG...\n");
        i64 total_edges_in_srg = 0;
        for (i32 i = 0; i < out.state_count - 1; i++) {
            for (i32 j = i + 1; j < out.state_count; j++) {
                // NOTE: We do not need to count the similarity of edges at this point, as we recount them later.
                // As this variant of the function is only called when testing, we accept this time waste.
                if (i > 0) { // Ensure that the SRG is connected.
                    i32 common_txs = count_common_txs(dfa, i, j);
                    if (common_txs < 128) continue;
                }

                count_edge(srg_edge_set, i, j);
                total_edges_in_srg++;
            }
        }

        finalize_adj_list(srg_edge_set);
    } else {
        // printf("  Constructing sparse SRG...\n");
        construct_sparse_srg_in_rounds(dfa, params, srg_edge_set, srg_time_info);
    }

    const clock_t t1 = clock();

    AdjacencyList tree_edge_set = make_adj_list(out.state_count);

    // Find the middle of an actual MST.
    // printf("  Constructing temporary MST.\n");
    {
        i64 tree_edge_count = 0;
        reset_adj_list(tree_edge_set);

        BitArray nodes_in_tree = ba_make(out.state_count);
        ba_zero(nodes_in_tree);

        std::priority_queue<PrimEdge, Vector<PrimEdge>, OrderPrimOnlyByWeight> queue;

        u32 prim_root = 0; // Pick arbitrary root.
        prim_add_edges_from_node(prim_root, 0, srg_edge_set, dfa, queue, nodes_in_tree);
        ba_set(nodes_in_tree, prim_root, true);

        // printf("   Running Prims.\n");
        i64 tree_weight = 0;
        i64 req_edges = out.state_count - 1;
        while (!queue.empty() && tree_edge_count < req_edges) {
            PrimEdge edge = queue.top();
            queue.pop();

            // Check if edge would cause a cycle.
            if (ba_get(nodes_in_tree, edge.dest)) continue;
            ba_set(nodes_in_tree, edge.dest, true);

            count_edge(tree_edge_set, edge.source, edge.dest);

            tree_edge_count++;
            tree_weight += edge.weight;

            // Queue further edges.
            prim_add_edges_from_node(edge.dest, edge.dist, srg_edge_set, dfa, queue, nodes_in_tree);
        }

        // printf("   Finalizing adjacency list.\n");
        finalize_adj_list(tree_edge_set);

        ba_delete(nodes_in_tree);
    }

    // Find middle of tree.
    // printf("  Finding middle of temporary MST.\n");
    i32 middle_of_mst = -1;
    {
        TreeBFSNode *queue = alloc_array<TreeBFSNode>(out.state_count);

        i32 last_node_hit_idx = bfs_to_exhaustion(queue, out.state_count, 0, tree_edge_set);

        i32 last_node_hit = queue[last_node_hit_idx].node;
        i32 head_idx = bfs_to_exhaustion(queue, out.state_count, last_node_hit, tree_edge_set);

        // Pick middle of path.
        // NOTE: This is a silly way to get the path... Could use two passes instead?
        Vector<i32> path;
        while (head_idx >= 0) {
            path.push_back(queue[head_idx].node);
            head_idx = queue[head_idx].parent;
        }

        middle_of_mst = path[path.size() / 2];

        free_and_reset_array(queue);
    }

    // Run modified Prims instead.
    i64 tree_edge_count = 0;
    {
        reset_adj_list(tree_edge_set);

        std::priority_queue<PrimEdge, Vector<PrimEdge>, OrderPrimHeuristic> queue;

        BitArray nodes_in_tree = ba_make(out.state_count);
        ba_zero(nodes_in_tree);

        // printf("   Running Prims.\n");

        // NOTE: The compression is sensitive to a good root. We don't yet know how to pick one.
        assert(middle_of_mst >= 0);
        u32 prim_root = (u32)middle_of_mst;
        prim_add_edges_from_node(prim_root, 0, srg_edge_set, dfa, queue, nodes_in_tree);
        ba_set(nodes_in_tree, prim_root, true);

        tree_edge_count = 0;
        i64 req_edges = out.state_count - 1;
        while (!queue.empty() && tree_edge_count < req_edges) {
            PrimEdge edge = queue.top();
            queue.pop();

            // Check if edge would cause a cycle.
            // We perform Find by iterating to the root of the tree.
            if (ba_get(nodes_in_tree, edge.dest)) continue;
            ba_set(nodes_in_tree, edge.dest, true);

            count_edge(tree_edge_set, edge.source, edge.dest);

            tree_edge_count++;

            // Queue further edges.
            prim_add_edges_from_node(edge.dest, edge.dist, srg_edge_set, dfa, queue, nodes_in_tree);
        }

        // printf("   Finalizing adjacency list.\n");
        finalize_adj_list(tree_edge_set);

        ba_delete(nodes_in_tree);
    }

    // Clean up as we go.
    free_ajd_list(srg_edge_set);

    // Construct default txs
    // printf("  Constructing default txs.\n");
    construct_default_txs_from_outside_going_in(out, tree_edge_set, tree_edge_count);

    clock_t t2 = clock();

    // Traverse nodes in order of increasing height and cut when a node would violate the diameter bound.
    if (diameter_bound > 0) {
        // Find root.
        // printf("  Finding root.\n");
        u32 tree_root = kNoState;
        for (i32 i = 0; i < out.state_count; i++) {
            if (out.default_txs[i] == kNoState) {
                assert(tree_root == kNoState);
                tree_root = i; // NOTE: We could break here but we don't to assert on double-assign.
            }
        }

        i32 *order  = alloc_array<i32>(out.state_count);
        i32 *height = alloc_array<i32>(out.state_count);
        i32 *tall   = alloc_array<i32>(out.state_count);

        memzero(height, out.state_count);

        // Construct order of increasing height.
        i32 order_index = 0;
        construct_preorder(tree_root, order_index, order, out.default_txs, tree_edge_set);
        assert(order_index == out.state_count);

        // printf("  Cutting default txs.\n");
        i64 cut_count  = 0;
        i64 cut_weight = 0;
        for (i32 i = 0; i < out.state_count; i++) { // Not including root.
            const u32 v = order[out.state_count - 1 - i]; // Smallest heights first.
            if (v == tree_root) continue;

            const u32 u = out.default_txs[v];
            assert(u != kNoState);

            i32 h = height[v] + 1;
            if (h + height[u] > diameter_bound) {
                if (h > height[u]) {
                    cut_count++;
                    cut_weight += count_common_txs(dfa, u, v); // NOTE: Slow to re-count just for stats...

                    out.default_txs[v] = kNoState;

                    delete_edge(tree_edge_set, u, v);
                    tree_edge_count--;
                } else {
                    const u32 t = tall[u];

                    cut_count++;
                    cut_weight += count_common_txs(dfa, u, t); // NOTE: Slow to re-count just for stats...

                    out.default_txs[t] = kNoState;

                    delete_edge(tree_edge_set, u, t);
                    tree_edge_count--;

                    height[u] = h;
                    tall[u]   = v;
                }
            } else if (h > height[u]) {
                height[u] = h;
                tall[u]   = v;
            }
        }

        // printf("  > Cut %d defaults txs of an avg weight %.2f\n", cut_count, (f32)cut_weight / (f32)cut_count);

        // NOTE: We need to re-orient the small trees.
        // printf("  Reconstructing default txs.\n");
        fill_array(out.default_txs, out.state_count, kNoState);
        construct_default_txs_from_outside_going_in(out, tree_edge_set, tree_edge_count);

        free_and_reset_array(tall);
        free_and_reset_array(height);
        free_and_reset_array(order);
    }

    free_ajd_list(tree_edge_set);

    const clock_t t3 = clock();

    copy_in_uncompressed_labelled_txs_with_counting(dfa, out);

    const clock_t t4 = clock();

    const f64 t_secs_total  = (f64)(t4 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_srg = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_mst = (f64)(t2 - t1) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_cut = (f64)(t3 - t2) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_txs = (f64)(t4 - t3) / (f64)CLOCKS_PER_SEC;

    const f64 t_secs_to_srg_hash  = (f64)srg_time_info.clocks_to_hash      / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_srg_edges = (f64)srg_time_info.clocks_to_add_edges / (f64)CLOCKS_PER_SEC;

    log(Log::Stat, "Time (secs): Total=%.4f\n\tSRGs=%.4f (%.2f%%)\n\t  Hash =%.4f (%.2f%%)\n\t  Edges=%.4f (%.2f%%)\n\tMST=%.4f (%.2f%%)\n\tCut=%.4f (%.2f%%)\n\tTxs =%.4f (%.2f%%)\n",
           t_secs_total,
           t_secs_to_srg,       100.0 * t_secs_to_srg / t_secs_total,
           t_secs_to_srg_hash,  100.0 * t_secs_to_srg_hash  / t_secs_total,
           t_secs_to_srg_edges, 100.0 * t_secs_to_srg_edges / t_secs_total,
           t_secs_to_mst,       100.0 * t_secs_to_mst / t_secs_total,
           t_secs_to_cut,       100.0 * t_secs_to_cut / t_secs_total,
           t_secs_to_txs,       100.0 * t_secs_to_txs / t_secs_total);
}

void calculate_node_depths(const DFA &dfa, i32 *depth) {
    fill_array(depth, dfa.state_count, -1);

    i32 *queue = alloc_array<i32>(dfa.state_count);
    i32 i0 = 0;
    i32 i1 = 1;

    queue[0] = 0; // Queue root.
    depth[0] = 0;
    while (i0 < dfa.state_count) {
        i32 cur = queue[i0++];
        for (i32 c = 0; c < dfa.alphabet_size; c++) {
            i32 dest = at(dfa.txs, cur, c);
            if (depth[dest] >= 0) continue;

            depth[dest] = depth[cur] + 1;
            queue[i1] = dest;
            i1++;
        }
    }
    assert(i1 == i0 && i0 == dfa.state_count);

    free_and_reset_array(queue);
}

void dfa_to_adfa_full_srg(const DFA &dfa, D2FA &out, const D2FAParams &param) {
    const i32 kMaxAlphabetSize = 1 << 8*sizeof(label_t);
    assert(dfa.alphabet_size <= kMaxAlphabetSize);

    out.alphabet_size = dfa.alphabet_size;
    out.state_count   = dfa.state_count;

    // NOTE: At this point we do not know the number of labelled transitions, so we do not allocate yet.
    // NOTE: We perform copies to construct a new D2FA without modifying the DFA.
    //       It might be more efficient to perform a transformation.
    out.tx_ranges        = alloc_array<Range>(out.state_count);
    out.default_txs      = alloc_array<state_t>(out.state_count);
    out.accepting_states = alloc_array<i32>(out.state_count);

    reset_mem_count();

    for (i32 i = 0; i < out.state_count; i++) out.default_txs[i] = kNoState;
    memcpy(out.accepting_states, dfa.accepting_states, sizeof(i32) * out.state_count);

    const clock_t t0 = clock();

    // Assign depth to all states in DFA by BFS.
    i32 *depth = alloc_array<i32>(out.state_count);
    calculate_node_depths(dfa, depth);

    // Get minimum spanning "tree" in the directed graph.
    const i32 adfa_k = 1; // Minimum number of depths to skip.
    i64 total_compressed_txs_count = 0;
    for (i32 i = 0; i < out.state_count; i++) {
        i32 max_dest   = -1;
        i32 max_weight = -1;
        for (i32 j = 0; j < out.state_count; j++) {
            if (depth[i] - depth[j] < adfa_k) continue;

            i32 common_txs = count_common_txs(dfa, i, j);
            if (common_txs > max_weight) {
                max_weight = common_txs;
                max_dest   = j;
            }
        }

        if (max_dest == -1) continue;
        assert(max_dest >= 0);
        out.default_txs[i] = max_dest;
        total_compressed_txs_count += max_weight;
    }

    const clock_t t1 = clock();

    copy_in_uncompressed_labelled_txs(dfa, out, total_compressed_txs_count);
   
    const clock_t t2 = clock();

    const f64 t_secs_total  = (f64)(t2 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_srg = (f64)(t1 - t0) / (f64)CLOCKS_PER_SEC;
    const f64 t_secs_to_txs = (f64)(t2 - t1) / (f64)CLOCKS_PER_SEC;

    log(Log::Stat, "Time (secs): Total=%.4f\n\tSRGs=%.4f (%.2f%%)\n\tTxs =%.4f (%.2f%%)\n",
           t_secs_total,
           t_secs_to_srg,    100.0 * t_secs_to_srg    / t_secs_total,
           t_secs_to_txs,    100.0 * t_secs_to_txs    / t_secs_total);

    free_and_reset_array(depth);
}


void dfa_to_adfa_sparse(const DFA &dfa, D2FA &out, const D2FAParams &params) {
    const i32 kMaxAlphabetSize = 1 << 8*sizeof(label_t);
    assert(dfa.alphabet_size <= kMaxAlphabetSize);

    out.alphabet_size = dfa.alphabet_size;
    out.state_count   = dfa.state_count;

    // NOTE: At this point we do not know the number of labelled transitions, so we do not allocate yet.
    // NOTE: We perform copies to construct a new D2FA without modifying the DFA.
    //       It might be more efficient to perform a transformation.
    out.tx_ranges        = alloc_array<Range>(out.state_count);
    out.default_txs      = alloc_array<state_t>(out.state_count);
    out.accepting_states = alloc_array<i32>(out.state_count);

    reset_mem_count();

    for (i32 i = 0; i < out.state_count; i++) out.default_txs[i] = kNoState;
    memcpy(out.accepting_states, dfa.accepting_states, sizeof(i32) * out.state_count);

    const i32 adfa_k = 1; // Minimum number of depths to skip.

    // Assign depth to all states in DFA by BFS.
    i32 *depth = alloc_array<i32>(out.state_count);
    calculate_node_depths(dfa, depth);

    i32 *saved = alloc_array<i32>(out.state_count);
    fill_array(saved, out.state_count, 0);

    // NOTE: Copy pasted from  paste.
    SRGTimingInfo timing_info;
    {
        u32 min_hash_k = params.min_hash_k;
        MinHashBoostApproach boost_approach = params.hash_approach;
        u32 min_hash_rounds = params.min_hash_r;

        i32 *hashes = alloc_array<i32>(dfa.state_count); // Hash for each state.

        Vector<i32> bucket_hashes;   // All hash keys of non-empty buckets.
        UnorderedMap<i32, Vector<state_t>> buckets; // Map from hash to bucket contents.

        for (u32 r = 0; r < min_hash_rounds; r++) {
            clock_t t0 = clock();

            if      (boost_approach == MinHashBoostApproach::KElements)         hash_states_k_smallest(dfa, min_hash_k, hashes);
            else if (boost_approach == MinHashBoostApproach::KPermutations)     hash_states_k_combined(dfa, min_hash_k, hashes);
            else if (boost_approach == MinHashBoostApproach::KSamplesReplace)   hash_states_k_samples_with_replace(dfa, min_hash_k, hashes);
            else if (boost_approach == MinHashBoostApproach::KSamplesNoReplace) hash_states_k_samples_no_replace(dfa, min_hash_k, hashes);
            else assert(false);

            clock_t t1 = clock();

            // Construct the buckets.
            bucket_hashes.clear();
            buckets.clear();
            for (i32 i = 0; i < dfa.state_count; i++) {
                const i32 h = hashes[i];

                Vector<state_t> &bucket = buckets[h];
                if (bucket.empty()) bucket_hashes.push_back(h);
                bucket.push_back(i);
            }

            // For each bucket, construct the edges.
            for (i32 b : bucket_hashes) {
                const Vector<state_t> &bucket = buckets[b];
                const u32 bucket_size = bucket.size();

                assert(bucket_size >= 1);
                if (bucket_size == 1) continue;

                for (u32 i = 0; i < bucket_size; i++) {
                    u32 j = (i + 1 + (random_u32() % (bucket_size - 1))) % bucket_size;
                    if (depth[bucket[i]] - depth[bucket[j]] < adfa_k) continue; // Wrong direction.
                    if (out.default_txs[bucket[i]] == bucket[j]) continue; // Same edge.

                    i32 common_txs = count_common_txs(dfa, bucket[i], bucket[j]);
                    if (common_txs <= saved[bucket[i]]) continue; // Saving nothing.

                    out.default_txs[bucket[i]] = bucket[j];
                    saved[bucket[i]] = common_txs;
                }
            }

            clock_t t2 = clock();
            timing_info.clocks_to_hash      += t1 - t0;
            timing_info.clocks_to_add_edges += t2 - t1;
        }

        free_and_reset_array(hashes);
    }

    copy_in_uncompressed_labelled_txs_with_counting(dfa, out);

    free_and_reset_array(depth);
    free_and_reset_array(saved);
}

typedef std::pair<state_t, state_t> StatePair;

struct PairHash
{
    size_t operator() (const std::pair<state_t, state_t> &pair) const {
        return ((pair.first + 60251L) * 11903L) + (pair.second * 100673L);
    }
};

struct TxLine {
    state_t u;
    state_t txs[256];
    bool accepting;
};

DFA merge(const DFA &a, const DFA &b) {
    assert(a.alphabet_size == b.alphabet_size);
    assert(a.alphabet_size == 256);

    DFA out;
    out.alphabet_size = a.alphabet_size;

    i64 k = 0;

    std::unordered_map<StatePair, state_t, PairHash> ids;
    ids.reserve(a.state_count + b.state_count);
    ids[{0, 0}] = k++;

    std::vector<TxLine> tx_lines;

    std::deque<StatePair> queue = { { 0, 0 } };
    while (!queue.empty()) {
        StatePair cur = queue.front();
        queue.pop_front();

        tx_lines.push_back({});
        TxLine &line = tx_lines.back();
        line.u = ids[cur];
        line.accepting =
            a.accepting_states[cur.first] > 0 ||
            b.accepting_states[cur.second] > 0;

        for (i32 i = 0; i < out.alphabet_size; i++) {
            StatePair dst = { at(a.txs, cur.first, i), at(b.txs, cur.second, i) };
            if (ids.count(dst) == 0) {
                queue.push_back(dst);
                ids[dst] = k++;
            }
            line.txs[i] = ids[dst];
        }
    }

    out.state_count = k;

    // Copy in transitions.
    out.accepting_states = alloc_array<i32>(out.state_count);
    out.txs = alloc_matrix<state_t>(a.alphabet_size, out.state_count);
    for (state_t i = 0; i < tx_lines.size(); i++) {
        assert(tx_lines[i].u == i);
        for (i32 j = 0; j < out.alphabet_size; j++) {
            at(out.txs, i, j) = tx_lines[i].txs[j];
        }
        out.accepting_states[i] = tx_lines[i].accepting;
    }

    return out;
}

state_t match_character(const DFA &f, state_t cur_state, label_t c) {
    return at(f.txs, (i32)cur_state, (i32)c);
}

// Match a zero-terminated string in the DFA, returning "true" iff a match is found.
bool match(const DFA &dfa, const char *str) {
    const char *c = str;
    state_t cur_state = 0; // Intial state.
    while (*c != 0) { // NOTE: We exit on zero character.
        cur_state = match_character(dfa, cur_state, *c++);
    }
    return dfa.accepting_states[cur_state] > 0;
}

state_t match_character(const D2FA &f, state_t cur_state, label_t c) {
    while (true) {
        const Range &tx_range = f.tx_ranges[cur_state];

        // Attempt to match a character.
        for (i32 i = 0; i < tx_range.count; i++) {
            const LabeledTx  &tx = f.txs[tx_range.index + i];
            if (tx.label == c) {
                return tx.dest;
            }
        }

        // Did not match and return, so follow default tx.
        assert(f.default_txs[cur_state] != kNoState); // Otherwise the D2FA is incomplete.
        cur_state = f.default_txs[cur_state];
    }

    return cur_state;
}

// Match a zero-terminated string in the D2FA, returning "true" iff a match is found.
bool match(const D2FA &f, const char *str) {
    const char *c = str;
    state_t cur_state = 0; // Intial state.
    while (*c != 0) { // NOTE: We exit on zero character.
        cur_state = match_character(f, cur_state, *c++);
    }
    return f.accepting_states[cur_state] > 0;
}

void free_dfa(DFA &dfa) {
    free_matrix(dfa.txs);
    free_and_reset_array(dfa.accepting_states);
}

void free_d2fa(D2FA &d2fa) {
    free_and_reset_array(d2fa.txs);
    free_and_reset_array(d2fa.tx_ranges);
    free_and_reset_array(d2fa.default_txs);
    free_and_reset_array(d2fa.accepting_states);
}
