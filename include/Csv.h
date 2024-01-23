#pragma once

#include <stddef.h>
#include <sys/types.h>

#include "Types.h"

/**
 * @brief Parse the height of a matrix stored as CSV
 *
 * @param[in] path Path to CSV file
 * @return Matrix height on success,
 *         <0 on error
 */
ssize_t csv_parse_matrix_height(const char* path); 

/**
 * @brief Parse the width of a matrix stored as CSV
 *
 * @param[in] path Path to CSV file
 * @return Matrix width on success,
 *         <0 on error
 */
ssize_t csv_parse_matrix_width(const char* path); 

/**
 * @brief Parse a vector stored as CSV
 *
 * @param[in] path Path to CSV file
 * @param[in] n Length of vector
 * @param[out] res Array with length n for storing vector
 * @return 0 on success,
 *         <0 on error
 */
int csv_parse_vector(const char* path, size_t n, real_t res[n]); 

/**
 * @brief Parse a matrix stored as CSV
 *
 * @param[in] path Path to CSV file
 * @param[in] m Height of matrix
 * @param[in] n Width of matrix
 * @param[out] res Flattened m by n array for storing matrix
 * @return 0 on success,
 *         <0 on error
 */
int csv_parse_matrix(const char* path, size_t m, size_t n, real_t res[m][n]); 

/**
 * @brief Store a vector as a CSV column
 *
 * @param[in] path Path for CSV file
 * @param[in] n Length of vector
 * @param[in] vec Array with length n to store
 * @return Number of bytes written if successful,
 *         <0 on error
 */
ssize_t csv_save_vector(const char* path, size_t n, const real_t vec[n]); 

/**
 * @brief Store a matrix as CSV
 *
 * @param[in] path Path for CSV file
 * @param[in] m Matrix height
 * @param[in] n Matrix width
 * @param[in] mat Flattened m by n array to store
 * @return Number of bytes written if successful,
 *         <0 on error
 */
ssize_t csv_save_matrix(const char* path, size_t m, size_t n, const real_t mat[m][n]); 
