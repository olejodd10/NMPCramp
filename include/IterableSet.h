#pragma once

#include <stdint.h>
#include <stddef.h>

// Struct definition/size is in header file to enable stack allocations.
typedef struct {
    size_t capacity; // Maximum number of elements to be stored in set
    size_t size; // Current size of set
    uint8_t* elements; // Set membership array
    size_t* next; // Next element in iteration sequence
    size_t* prev; // Previous element in iteration sequence
    size_t first; // First element in iteration sequence
    size_t last; // Last element in iteration sequence
} iterable_set_t;

/**
 * @brief Initialize iterable set.
 *
 * @param[out] set Iterable set instance
 * @param[in] capacity Maximum number of elements to be stored in set
 */
void iterable_set_init(iterable_set_t* set, size_t capacity);

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
 * @brief Checks whether element is in iterable set.
 *
 * @param[in] set Iterable set instance
 * @param[in] element Element to check for
 * @return 1 if element is in set, 0 otherwise
 */
uint8_t iterable_set_contains(const iterable_set_t* set, size_t element);

/**
 * @brief Returns the first element of a set iteration.
 *
 * @param[in] set Iterable set instance
 * @return The first element of a set iteration
 */
size_t iterable_set_first(const iterable_set_t* set);

/**
 * @brief Returns the next element in a set iteration.
 *
 * @param[in] set Iterable set instance
 * @param[in] element The previous element
 * @return The next element
 */
size_t iterable_set_next(const iterable_set_t* set, size_t element);

/**
 * @brief Returns the value returned by iterable_set_next when called with the last element.
 *
 * @param[in] set Iterable set instance
 * @return The value returned by iterable_set_next when called with the last element
 */
size_t iterable_set_end(const iterable_set_t* set);