#include "IndexedVectors.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "Types.h"

void indexed_vectors_init(indexed_vectors_t *vecs, size_t capacity, size_t length, size_t max_index) {
    vecs->capacity = capacity;
    vecs->length = length;
    vecs->max_index = max_index;
	vecs->where_is = (size_t*)malloc(max_index*sizeof(size_t)); 
	vecs->whats_here = (size_t*)malloc(capacity*sizeof(size_t)); 
    vecs->values = (real_t*)malloc(capacity*length*sizeof(real_t));
    for (size_t i = 0; i < vecs->max_index; ++i) {
        vecs->where_is[i] = vecs->capacity; 
    }
    for (size_t i = 0; i < vecs->capacity; ++i) {
        vecs->whats_here[i] = vecs->max_index; 
    }
    vecs->size = 0;
}

void indexed_vectors_destroy(indexed_vectors_t *vecs) {
    free(vecs->where_is);
    free(vecs->whats_here);
    free(vecs->values);
}

void indexed_vectors_clear(indexed_vectors_t *vecs) {
    for (size_t i = 0; i < vecs->capacity; ++i) {
        if (vecs->whats_here[i] != vecs->max_index) {
            vecs->where_is[vecs->whats_here[i]] = vecs->capacity;
            vecs->whats_here[i] = vecs->max_index; 
        }
    }
    vecs->size = 0;
}

static size_t first_available_slot(const indexed_vectors_t *vecs) {
    for (size_t i = 0; i < vecs->capacity; ++i) {
        if (vecs->whats_here[i] == vecs->max_index) {
            return i;
        }
    }
    return vecs->capacity;
}

real_t* indexed_vectors_get_mut(indexed_vectors_t *vecs, size_t index) {
    if (index >= vecs->max_index || vecs->where_is[index] == vecs->capacity) {
        return NULL;
    } else {
        return &vecs->values[vecs->length*vecs->where_is[index]];
    }
}

const real_t* indexed_vectors_get(const indexed_vectors_t *vecs, size_t index) {
    if (index >= vecs->max_index || vecs->where_is[index] == vecs->capacity) {
        return NULL;
    } else {
        return &vecs->values[vecs->length*vecs->where_is[index]];
    }
}

void indexed_vectors_insert(indexed_vectors_t *vecs, size_t index, const real_t *vec) {
    if (index >= vecs->max_index || vecs->where_is[index] != vecs->capacity 
            || vecs->size == vecs->capacity) {
        return;
    }
    size_t target = first_available_slot(vecs);
    vecs->whats_here[target] = index;
    vecs->where_is[index] = target;
    memcpy(&vecs->values[vecs->length*target], vec, vecs->length*sizeof(real_t));
    vecs->size++;
}

void indexed_vectors_remove(indexed_vectors_t *vecs, size_t index) {
    if (index >= vecs->max_index || vecs->where_is[index] == vecs->capacity) {
        return;
    }
    vecs->whats_here[vecs->where_is[index]] = vecs->max_index;
    vecs->where_is[index] = vecs->capacity;
    vecs->size--;
}
