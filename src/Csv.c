#include "Csv.h"

#include <stddef.h>
#include <sys/types.h>
#include <stdio.h>

#include "Types.h"

#define SAVE_FORMAT "%.4le"

static size_t parse_matrix_height(FILE* f) {
    size_t count = 0;
    fpos_t pos;
    fgetpos(f, &pos);
    rewind(f);
    while (fscanf(f, "%*s") != EOF) {
        ++count;
    }
    fsetpos(f, &pos);
    return count;
}

static size_t parse_matrix_width(FILE* f) {
    size_t count = 0;
    fpos_t pos;
    fgetpos(f, &pos);
    rewind(f);
    while (fscanf(f, "%*f,") != EOF) {
        ++count;
        int c = fgetc(f);
        if (c == '\n' || c =='\r') {
            break;
        }
        ungetc(c, f);
    }
    fsetpos(f, &pos);
    return count;
}

static int assert_size(FILE* f, size_t m, size_t n) {
    if (m != parse_matrix_height(f) || n != parse_matrix_width(f)) {
        return -1;
    }
    return 0;
}

// Accepts both row and column vectors
static int parse_row(FILE* f, size_t n, real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        int ret = fscanf(f, "%f,", &res[i]);
        if (ret != 1) {
            return -1;
        }
    }
    return 0;
}

ssize_t csv_parse_matrix_height(const char* path) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    size_t ret = parse_matrix_height(f);
    fclose(f);
    return ret;
}

ssize_t csv_parse_matrix_width(const char* path) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    size_t ret = parse_matrix_width(f);
    fclose(f);
    return ret;
}


int csv_parse_vector(const char* path, size_t n, real_t res[n]) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    int ret = assert_size(f, n, 1); // CSVs with one row could just as well be accepted
    if (ret < 0) {
        return ret;
    }
    ret = parse_row(f, n, res);
    fclose(f);
    return ret;
}

int csv_parse_matrix(const char* path, size_t m, size_t n, real_t res[m][n]) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    int ret = assert_size(f, m, n);
    if (ret < 0) {
        return ret;
    }
    for (size_t i = 0; i < m; ++i) {
        ret = parse_row(f, n, res[i]);
        if (ret < 0) {
            return ret;
        }
    }
    fclose(f);
    return ret;
}

static ssize_t write_row(FILE* f, size_t n, const real_t vec[n]) {
    size_t sum = 0;
    int ret = 0;
    for (size_t i = 0; i < n-1; ++i) {
        ret = fprintf(f, SAVE_FORMAT ",", vec[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = fprintf(f, SAVE_FORMAT "\n", vec[n-1]);
    if (ret < 0) {
        return ret;
    } else {
        sum += ret;
    }
    return sum;
}

// Writes as a column vector
ssize_t csv_save_vector(const char* path, size_t n, const real_t vec[n]) {
    FILE* f = fopen(path, "w+");
    if (f == NULL) {
        return -1;
    }
    size_t sum = 0;
    int ret;
    for (size_t i = 0; i < n; ++i) {
        ret = fprintf(f, SAVE_FORMAT "\n", vec[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = assert_size(f, n, 1);
    if (ret < 0) {
        return ret;
    }
    fclose(f);
    return sum;
}

ssize_t csv_save_matrix(const char* path, size_t m, size_t n, const real_t mat[m][n]) {
    FILE* f = fopen(path, "w+");
    if (f == NULL) {
        return -1;
    }
    size_t sum = 0;
    int ret;
    for (size_t i = 0; i < m; ++i) {
        ret = write_row(f, n, mat[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = assert_size(f, m, n);
    if (ret < 0) {
        return ret;
    }
    fclose(f);
    return sum;
}
