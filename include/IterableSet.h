#pragma once

#include <stdint.h>
#include <stddef.h>

// Struct definition/size is in header file to enable stack allocations.
typedef struct {
    size_t capacity; // Maximum number of elements to be stored in set
    size_t size; // Current size of set
    size_t pivot; // Index of first element that is counted as part of partition if in set
    size_t partition; // Number of elements >=pivot in set
    uint8_t* elements; // Set membership array
    size_t* ordering; // Elements in given positions
    size_t* ordering_of; // Positions of given elements
} iterable_set_t;

/**
 * @brief Initialize iterable set.
 *
 * @param[out] set Iterable set instance
 * @param[in] capacity Maximum number of elements to be stored in set
 * @param[in] pivot Index of first element that is counted as part of partition if in set
 */
void iterable_set_init(iterable_set_t* set, size_t capacity, size_t pivot);

/**
 * @brief Destroy iterable set.
 *
 * @param[out] set Iterable set instance
 */
void iterable_set_destroy(iterable_set_t* set);

/**
 * @brief Remove all elements from iterable set.
 *
 * @param[out] set Iterable set instance
 */
void iterable_set_clear(iterable_set_t* set);

/**
 * @brief Insert element into iterable set.
 *
 * @param[out] set Iterable set instance
 * @param[in] element Element to insert
 * @note does nothing if element is larger than set capacity or already in set
 */
void iterable_set_insert(iterable_set_t* set, size_t element);

/**
 * @brief Remove element from iterable set.
 *
 * @param[out] set Iterable set instance
 * @param[in] element Element to remove
 * @note does nothing if element is larger than set capacity or not in set
 */
void iterable_set_remove(iterable_set_t* set, size_t element);

/**
 * @brief Returns size of iterable set.
 *
 * @param[in] set Iterable set instance
 * @return Number of elements currently in set
 */
size_t iterable_set_size(const iterable_set_t* set);

/**
 * @brief Returns size of iterable set partition.
 *
 * @param[in] set Iterable set instance
 * @return Number of elements currently in set partition
 */
size_t iterable_set_partition(const iterable_set_t* set);

/**
 * @brief Checks whether element is in iterable set.
 *
 * @param[in] set Iterable set instance
 * @param[in] element Element to check for
 * @return 1 if element is in set, 0 otherwise
 */
uint8_t iterable_set_contains(const iterable_set_t* set, size_t element);

size_t iterable_set_nth(const iterable_set_t *set, size_t n);

size_t iterable_set_nth_low(const iterable_set_t *set, size_t n);

size_t iterable_set_nth_high(const iterable_set_t *set, size_t n);

size_t iterable_set_nth_missing_low(const iterable_set_t *set, size_t n);

size_t iterable_set_nth_missing_high(const iterable_set_t *set, size_t n);

/**
 * @brief Returns the value returned by iterable_set_nth* when called with an overshooting index.
 *
 * @param[in] set Iterable set instance
 * @return The value returned by iterable_set_nth* when called with an invalid index
 */
size_t iterable_set_end(const iterable_set_t* set);
