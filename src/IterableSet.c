#include "IterableSet.h"

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

void iterable_set_init(iterable_set_t* set, size_t capacity) {
    set->capacity = capacity;
    set->elements = (uint8_t*)malloc(capacity*sizeof(uint8_t));
    set->next = (size_t*)malloc(capacity*sizeof(size_t));
    set->prev = (size_t*)malloc(capacity*sizeof(size_t));
    memset(set->elements, 0, capacity);
    set->first = capacity;
    set->last = capacity;
    set->size = 0;
}

void iterable_set_destroy(iterable_set_t* set) {
	free(set->elements);
	free(set->next);
	free(set->prev);
}

void iterable_set_clear(iterable_set_t* set) {
    for (size_t i = iterable_set_first(set); i != iterable_set_end(set); i = iterable_set_next(set,i)) {
        set->elements[i] = 0; // Simplified removal
    }
    set->first = set->capacity;
    set->last = set->capacity;
    set->size = 0;
}

void iterable_set_insert(iterable_set_t* set, size_t element) {
    if (element >= set->capacity || set->elements[element]) {
        return;
    }
    set->elements[element] = 1;
    if (set->last == set->capacity) { // First element
        set->first = element;
        set->prev[element] = set->capacity; // Can drop this unless iterating backwards
    } else {
        set->next[set->last] = element;
        set->prev[element] = set->last;
    }
    set->next[element] = set->capacity;
    set->last = element;
    set->size++;
}

void iterable_set_remove(iterable_set_t* set, size_t element) {
    if (element >= set->capacity || !set->elements[element]) {
        return;
    }
    // Don't need to update prev and next for element because they only need to be correct for elements currently in the set
    set->elements[element] = 0;
    if (set->first == element) {
        set->first = set->next[element];
        if (set->last == element) {
            set->last = set->capacity;
        } else {
            set->prev[set->next[element]] = set->capacity; // Can drop this unless iterating backwards
        }
    } else if (set->last == element) {
        set->last = set->prev[element];
        set->next[set->prev[element]] = set->capacity;
    } else {
        set->next[set->prev[element]] = set->next[element];
        set->prev[set->next[element]] = set->prev[element];
    }
    set->size--;
}

size_t iterable_set_size(const iterable_set_t* set) {
    return set->size;
}

uint8_t iterable_set_contains(const iterable_set_t* set, size_t element) {
    return element < set->capacity && set->elements[element];
}

size_t iterable_set_first(const iterable_set_t* set) {
    return set->first;
}

// Iterates through active elements in order of insertion
size_t iterable_set_next(const iterable_set_t* set, size_t element) {
    return element < set->capacity ? set->next[element] : set->capacity;
}

size_t iterable_set_end(const iterable_set_t* set) {
    return set->capacity;
}
