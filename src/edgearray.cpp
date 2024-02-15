#include "edgearray.h"

#include "util.h"

EdgeArray earray_init(i32 size) {
    EdgeArray a;
    a.count = 0;
    a.size  = size;
    a.data  = alloc_array<Edge>(size);
    return a;
}

void earray_free(EdgeArray &a) {
    free_and_reset_array(a.data);
    a.count = 0;
    a.size  = 0;
}

void earray_clear(EdgeArray &a) {
    a.count = 0;
}

void earray_push(EdgeArray &a, const Edge &e) {
    assert(a.count < a.size);
    a.data[a.count++] = e;
}
