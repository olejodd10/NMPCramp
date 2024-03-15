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
        size_t swap_element = set->ordering[set->pivot + set->partition];
        swap(set->ordering, set->pivot + set->partition, set->ordering_of[element]);
        swap(set->ordering_of, swap_element, element);
        set->partition++;
    } else {
        size_t swap_element = set->ordering[set->size - set->partition];
        swap(set->ordering, set->size - set->partition, set->ordering_of[element]);
        swap(set->ordering_of, swap_element, element);
    }
    set->elements[element] = 1;
    set->size++;
}

void iterable_set_remove(iterable_set_t* set, size_t element) {
    if (element >= set->capacity || !set->elements[element]) {
        return;
    }
    if (element >= set->pivot) {
        size_t swap_element = set->ordering[set->pivot + set->partition - 1];
        swap(set->ordering, set->pivot + set->partition - 1, set->ordering_of[element]);
        swap(set->ordering_of, swap_element, element);
        set->partition--;
    } else {
        size_t swap_element = set->ordering[set->size - set->partition - 1];
        swap(set->ordering, set->size - set->partition - 1, set->ordering_of[element]);
        swap(set->ordering_of, swap_element, element);
    }
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
    size_t num_low = set->size - set->partition;
    if (n < num_low) {
        return iterable_set_nth_low(set, n);
    } else {
        return iterable_set_nth_high(set, n - num_low);
    }
}

size_t iterable_set_nth_low(const iterable_set_t *set, size_t n) {
    if (n >= set->size - set->partition) {
        return set->capacity;
    } else {
        return set->ordering[n];
    }
}

size_t iterable_set_nth_high(const iterable_set_t *set, size_t n) {
    if (n >= set->partition) {
        return set->capacity;
    } else {
        return set->ordering[set->pivot + n];
    }
}

size_t iterable_set_nth_missing_low(const iterable_set_t *set, size_t n) {
    size_t num_low = set->size - set->partition;
    size_t num_missing_low = set->pivot - num_low;
    if (n >= num_missing_low) {
        return set->capacity;
    } else {
        return set->ordering[num_low + n];
    }
}

size_t iterable_set_nth_missing_high(const iterable_set_t *set, size_t n) {
    size_t num_high = set->partition;
    size_t num_missing_high = set->capacity - set->pivot - num_high;
    if (n >= num_missing_high) {
        return set->capacity;
    } else {
        return set->ordering[set->pivot + num_high + n];
    }
}

size_t iterable_set_end(const iterable_set_t* set) {
    return set->capacity;
}
