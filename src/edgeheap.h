#ifndef __EDGEHEAP_H_
#define __EDGEHEAP_H_

#include <vector>
#include <unordered_map>

#include "common.h"
#include "cntalloc.h"

struct HeapEdge {
    SmallEdge edge;
    i32  weight;
};

struct EdgeHeap {
    Vector<HeapEdge> edges;
    UnorderedMap<u64, i32> positions; // Key -> Position.
    i32 resize_count = 0;
};

EdgeHeap eheap_init(i64 size);
void eheap_free(EdgeHeap &heap);

void eheap_clear(EdgeHeap &heap);

void eheap_push(EdgeHeap &heap, SmallEdge &edge, i32 weight);
SmallEdge eheap_pop(EdgeHeap &heap);

void eheap_update_weight(EdgeHeap &heap, u64 edge_key, i32 new_weight);

inline bool eheap_is_empty(const EdgeHeap &heap) {
    return heap.edges.empty();
}
inline bool eheap_contains(EdgeHeap &heap, u64 edge_key) {
    return heap.positions.count(edge_key);
}

#endif // __EDGEHEAP_H_
