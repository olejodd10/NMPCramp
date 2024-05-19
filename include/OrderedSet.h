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
} ordered_set_t;

/**
 * @brief Initialize ordered set.
 *
 * @param[out] set Ordered set instance
 * @param[in] capacity Maximum number of elements to be stored in set
 * @param[in] pivot Index of first element that is counted as part of partition if in set
 */
void ordered_set_init(ordered_set_t* set, size_t capacity, size_t pivot);

/**
 * @brief Destroy ordered set.
 *
 * @param[out] set Ordered set instance
 */
void ordered_set_destroy(ordered_set_t* set);

/**
 * @brief Remove all elements from ordered set.
 *
 * @param[out] set Ordered set instance
 */
void ordered_set_clear(ordered_set_t* set);

/**
 * @brief Insert element into ordered set.
 *
 * @param[out] set Ordered set instance
 * @param[in] element Element to insert
 * @note does nothing if element is larger than set capacity or already in set
 */
void ordered_set_insert(ordered_set_t* set, size_t element);

/**
 * @brief Remove element from ordered set.
 *
 * @param[out] set Ordered set instance
 * @param[in] element Element to remove
 * @note does nothing if element is larger than set capacity or not in set
 */
void ordered_set_remove(ordered_set_t* set, size_t element);

/**
 * @brief Returns size of ordered set.
 *
 * @param[in] set Ordered set instance
 * @return Number of elements currently in set
 */
size_t ordered_set_size(const ordered_set_t* set);

/**
 * @brief Returns size of ordered set partition.
 *
 * @param[in] set Ordered set instance
 * @return Number of elements currently in set partition
 */
size_t ordered_set_partition(const ordered_set_t* set);

/**
 * @brief Checks whether element is in ordered set.
 *
 * @param[in] set Ordered set instance
 * @param[in] element Element to check for
 * @return 1 if element is in set, 0 otherwise
 */
uint8_t ordered_set_contains(const ordered_set_t* set, size_t element);

size_t ordered_set_nth(const ordered_set_t *set, size_t n);

/**
 * @brief Returns the value returned by ordered_set_nth when called with an overshooting index.
 *
 * @param[in] set Ordered set instance
 * @return The value returned by ordered_set_nth when called with an overshooting index
 */
size_t ordered_set_end(const ordered_set_t* set);
