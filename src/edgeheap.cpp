#include "edgeheap.h"

#include "util.h"

inline i32 parent(i32 index) {
    return (index - 1) / 2;
}

EdgeHeap eheap_init(i64 size) {
    EdgeHeap heap;
    heap.edges.reserve(size);
    heap.positions.reserve(size);
    return heap;
}

void eheap_free(EdgeHeap &heap) {
    // Does nothing due to STL use..
}

void eheap_clear(EdgeHeap &heap) {
    heap.edges.clear();
    heap.positions.clear();

    heap.resize_count = 0;
}

void swap(EdgeHeap &heap, i32 pos_a, i32 pos_b) {
    const SmallEdge &edge_a = heap.edges[pos_a].edge;
    const SmallEdge &edge_b = heap.edges[pos_b].edge;
    u64 key_a = to_edge_key(edge_a.source, edge_a.dest);
    u64 key_b = to_edge_key(edge_b.source, edge_b.dest);

    HeapEdge a = heap.edges[pos_a];
    heap.edges[pos_a] = heap.edges[pos_b];
    heap.edges[pos_b] = a;

    heap.positions[key_a] = pos_b;
    heap.positions[key_b] = pos_a;
}

void sift_up(EdgeHeap &heap, i32 pos) {
    i32 i = pos;
    while (i > 0 && heap.edges[i].weight > heap.edges[parent(i)].weight) {
        swap(heap, i, parent(i));
        i = parent(i);
    }
}

void sift_down(EdgeHeap &heap, i32 pos) {
    u32 heap_count = heap.edges.size();
    u32 i = pos;
    while (true) {
        u32 a = i*2 + 1;
        u32 b = i*2 + 2;
        if (a >= heap_count) break;

        bool choose_a = b >= heap_count || heap.edges[a].weight > heap.edges[b].weight;
        u32 m = choose_a ? a : b;
        bool should_move = heap.edges[m].weight > heap.edges[i].weight;
        if (!should_move) break;

        swap(heap, i, m);
        i = m;
    }
}

void eheap_push(EdgeHeap &heap, SmallEdge &edge, i32 weight) {
    i32 position = heap.edges.size();
    u64 edge_key = to_edge_key(edge.source, edge.dest);
    heap.positions[edge_key] = position;

    if (heap.edges.size() == heap.edges.capacity()) heap.resize_count++;
    heap.edges.push_back({ edge, weight });

    sift_up(heap, position);
}

SmallEdge eheap_pop(EdgeHeap &heap) {
    SmallEdge root = heap.edges[0].edge;

    u64 root_key = to_edge_key(root.source, root.dest);
    assert(heap.positions.count(root_key)); // NOTE: Sanity check, remove for speed.
    assert(heap.positions[root_key] == 0);

    swap(heap, 0, heap.edges.size() - 1);
    assert(heap.edges[heap.edges.size() - 1].edge.source == root.source);
    assert(heap.edges[heap.edges.size() - 1].edge.dest   == root.dest);
    heap.edges.pop_back();
    heap.positions.erase(root_key);

    sift_down(heap, 0);

    return root;
}

void eheap_update_weight(EdgeHeap &heap, u64 edge_key, i32 new_weight) {
    if (heap.positions.count(edge_key) == 0) return; // Only if we have the edge.

    i32 pos = heap.positions[edge_key];
    bool old_weight        = heap.edges[pos].weight;
    heap.edges[pos].weight = new_weight;

    if (new_weight > old_weight) {
        sift_up(heap, pos);
    } else if (new_weight < old_weight) {
        sift_down(heap, pos);
    }
}
