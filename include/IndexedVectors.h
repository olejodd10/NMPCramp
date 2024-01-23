#pragma once

#include <stddef.h>

#include "Types.h"

// Struct definition/size is in header file to enable stack allocations.
typedef struct {
    size_t length; // Length of vectors to be stored
    size_t capacity; // Number of slots for storing vectors
    size_t max_index; // Maximum index of vectors to be stored
    size_t size; // Current number of slots used
    size_t* whats_here; // The indices of vectors stored in given slots
    size_t* where_is; // The slots where vectors with given indices are stored
    real_t* values; // The stored vector data
} indexed_vectors_t;

/**
 * @brief Initialize indexed vectors instance.
 *
 * @param[out] vecs Indexed vectors instance
 * @param[in] capacity Maximum number of vectors to be stored in vecs
 * @param[in] length Length of vectors to be store in vecs
 * @param[in] max_index Maximum index of vectors to be stored in vecs
 */
void indexed_vectors_init(indexed_vectors_t *vecs, size_t capacity, size_t length, size_t max_index);

/**
 * @brief Destroy indexed vectors instance.
 *
 * @param[out] vecs Indexed vectors instance
 */
void indexed_vectors_destroy(indexed_vectors_t *vecs);

/**
 * @brief Remove all vectors from indexed vectors instance.
 *
 * @param[out] vecs Indexed vectors instance
 */
void indexed_vectors_clear(indexed_vectors_t *vecs);

/**
 * @brief Returns a pointer to a non const stored vector.
 *
 * @param[in] vecs Indexed vectors instance
 * @param[in] index Index of wanted vector
 * @return Pointer to non const vector if index is stored, NULL otherwise
 */
real_t* indexed_vectors_get_mut(indexed_vectors_t *vecs, size_t index);

/**
 * @brief Returns a pointer to a const stored vector.
 *
 * @param[in] vecs Indexed vectors instance
 * @param[in] index Index of wanted vector
 * @return Pointer to const vector if index is stored, NULL otherwise
 */
const real_t* indexed_vectors_get(const indexed_vectors_t *vecs, size_t index);

/**
 * @brief Store an indexed vector.
 *
 * @param[out] vecs Indexed vectors instance
 * @param[in] index Index of vector to insert
 * @param[in] vec Vector to insert
 * @note Does nothing if index exceeds max_index, if index is already stored or if capacity is met
 * @warning vec must have length equal to the parameter length passed to indexed_vectors_init
 */
void indexed_vectors_insert(indexed_vectors_t *vecs, size_t index, const real_t *vec);

/**
 * @brief Remove an indexed vector.
 *
 * @param[out] vecs Indexed vectors instance
 * @param[in] index Index of vector to remove
 * @note Does nothing if index is not already stored
 */
void indexed_vectors_remove(indexed_vectors_t *vecs, size_t index);
