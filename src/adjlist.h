#ifndef ADJ_LIST_H_
#define ADJ_LIST_H_

#include <vector>

#include "common.h"
#include "util.h"

const i32 kDeletedMarker = -1;

//  NOTE: This is undirected!
struct SetEdge {
    i32 a;
    i32 b;
};

// NOTE: This is undirected!
struct AdjacencyList {
    i64 node_count; // Number of nodes in the graph.

    i64 *starts; // [i] = Index of first neighbour (in dests) for node i.
    i64 *dests;
    Vector<SetEdge> buffer; // All edge destinations, with all neighbour of each node being contiguous.
};

AdjacencyList make_adj_list(i32 node_count);

void reset_adj_list(AdjacencyList &set);
void finalize_adj_list(AdjacencyList &set);

void free_ajd_list(AdjacencyList &set);

void delete_neighbour_one_way(AdjacencyList &set, i32 node, i32 target);

inline void count_edge(AdjacencyList &set, i32 a, i32 b) {
    set.starts[a]++;
    set.starts[b]++;
    set.buffer.push_back(SetEdge { a, b });
}

inline void get_edge_range(const AdjacencyList &edge_set, i32 node, i64 &out_a, i64 &out_b) {
    out_a = edge_set.starts[node];
    out_b = edge_set.starts[node + 1];
    while (edge_set.dests[out_b - 1] == kDeletedMarker) out_b--;
}

inline void delete_edge(AdjacencyList &set, i32 a, i32 b) {
    delete_neighbour_one_way(set, a, b);
    delete_neighbour_one_way(set, b, a);
}

#endif /* ADJ_LIST_H_ */
