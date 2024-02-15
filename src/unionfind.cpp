#include "unionfind.h"

#include "util.h"

UnionFind uf_init(i32 count) {
    UnionFind uf;
    uf.count = count;
    uf.szs = alloc_array<i32>(count);
    uf.ids = alloc_array<i32>(count);

    uf_reset(uf);
    return uf;
}

void uf_delete(UnionFind &uf) {
    free_and_reset_array(uf.ids);
    free_and_reset_array(uf.szs);
}

void uf_reset(UnionFind &uf) {
   for (i32 i = 0; i < uf.count; i++) {
        uf.ids[i] = i;
        uf.szs[i] = 1;
    }
}

i32 uf_find(UnionFind &uf, i32 a) {
    while (a != uf.ids[a]) {
        uf.ids[a] = uf.ids[uf.ids[a]]; // Path compression.
        a = uf.ids[a];
    }
    return a;
}

void uf_union(UnionFind &uf, i32 a, i32 b) {
    a = uf_find(uf, a);
    b = uf_find(uf, b);
    if (uf.szs[a] < uf.szs[b]) {
        uf.ids[a] = b;
        uf.szs[b] += uf.szs[a];
    } else {
        uf.ids[b]  = a;
        uf.szs[a] += uf.szs[b];
    }
}
