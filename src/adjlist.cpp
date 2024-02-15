#include "adjlist.h"

#include "util.h"

AdjacencyList make_adj_list(i32 node_count) {
    AdjacencyList set;
    set.node_count = node_count;

    set.starts = alloc_array<i64>(node_count + 1); // Space for "beyond-the-end" element.
    set.dests  = NULL;

    memzero(set.starts, set.node_count + 1);

    set.buffer.reserve(2 * node_count); // Minimum for a connected graph.
    return set;
}

void free_ajd_list(AdjacencyList &set) {
    free_and_reset_array(set.starts);
    free_and_reset_array(set.dests);
}

void reset_adj_list(AdjacencyList &set) {
    memzero(set.starts, set.node_count + 1);
    free_and_reset_array(set.dests);
    set.buffer.clear();
}

void finalize_adj_list(AdjacencyList &set) {
    const i32 starts_size = set.node_count + 1;

    // We calculate one final "start index", at the very end of the array, which is actually the end of the final range.
    // This simplifies the logic when we want to know the end (and size) of any interval, as we can always just look at the next "start".
    // To actually fill the ranges we use an additional array of pointers, as we need to keep the start indices for afterwards.
    i64 *ptrs = alloc_array<i64>(starts_size);

    prefix_sum<i64>(set.starts, starts_size);
    memcpy(ptrs, set.starts, sizeof(i64) * starts_size);
    assert(set.starts[starts_size - 1] == 2 * (i64)set.buffer.size()); // Should have all edges twice.

    // Construct the neighbour lists.
    assert(set.dests == NULL); // Should only be filled once - now.
    set.dests = alloc_array<i64>(2 * set.buffer.size());
    for (const SetEdge &edge : set.buffer) {
        set.dests[ptrs[edge.a]] = edge.b;
        set.dests[ptrs[edge.b]] = edge.a;
        ptrs[edge.a]++;
        ptrs[edge.b]++;
    }
    set.buffer.clear(); // NOTE: We don't need to do this, but I think it reduces confusion.

    free_and_reset_array(ptrs);
}

void delete_neighbour_one_way(AdjacencyList &set, i32 node, i32 target) {
    i64 l, r;
    get_edge_range(set, node, l, r);
    for (i32 i = l; i < r; i++) {
        if (set.dests[i] == target) {
            swap_in_array(&set.dests[l], r - l - 1, i - l);
            set.dests[r - 1] = kDeletedMarker;
            return;
        }
    }
    assert(false);
}
