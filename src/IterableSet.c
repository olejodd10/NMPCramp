#include "IterableSet.h"

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

static void swap(size_t *arr, size_t i, size_t j) {
    size_t temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

void iterable_set_init(iterable_set_t* set, size_t capacity, size_t pivot) {
    set->capacity = capacity;
    set->pivot = pivot;
    set->elements = (uint8_t*)malloc(capacity*sizeof(uint8_t));
    set->ordering = (size_t*)malloc(capacity*sizeof(size_t));
    set->ordering_of = (size_t*)malloc(capacity*sizeof(size_t));
    memset(set->elements, 0, capacity);
    for (size_t i = 0; i < capacity; ++i) {
        set->ordering[i] = i;
        set->ordering_of[i] = i;
    }
    set->size = 0;
    set->partition = 0;
}

void iterable_set_destroy(iterable_set_t* set) {
	free(set->elements);
	free(set->ordering);
	free(set->ordering_of);
}

void iterable_set_clear(iterable_set_t* set) {
    memset(set->elements, 0, set->capacity);
    for (size_t i = 0; i < set->capacity; ++i) {
        set->ordering[i] = i;
        set->ordering_of[i] = i;
    }
    set->size = 0;
    set->partition = 0;
}

void iterable_set_insert(iterable_set_t* set, size_t element) {
    if (element >= set->capacity || set->elements[element]) {
        return;
    }
    if (element >= set->pivot) {
        set->partition++;
    } 
    size_t swap_element = set->ordering[set->size];
    swap(set->ordering, set->size, set->ordering_of[element]);
    swap(set->ordering_of, swap_element, element);
    set->elements[element] = 1;
    set->size++;
}

void iterable_set_remove(iterable_set_t* set, size_t element) {
    if (element >= set->capacity || !set->elements[element]) {
        return;
    }
    if (element >= set->pivot) {
        set->partition--;
    } 
    size_t swap_element = set->ordering[set->size - 1];
    swap(set->ordering, set->size - 1, set->ordering_of[element]);
    swap(set->ordering_of, swap_element, element);
    set->elements[element] = 0;
    set->size--;
}

size_t iterable_set_size(const iterable_set_t* set) {
    return set->size;
}

size_t iterable_set_partition(const iterable_set_t* set) {
    return set->partition;
}

uint8_t iterable_set_contains(const iterable_set_t* set, size_t element) {
    return element < set->capacity && set->elements[element];
}

size_t iterable_set_nth(const iterable_set_t *set, size_t n) {
    if (n >= set->size) {
        return set->capacity;
    } else {
        return set->ordering[n];
    }
}

size_t iterable_set_end(const iterable_set_t* set) {
    return set->capacity;
}
