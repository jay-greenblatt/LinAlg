#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, double low, double high) {
    srand(42);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows <= 0 || cols <= 0) {
      PyErr_SetString(PyExc_TypeError, "Invalid Dimensions!");
      return -1;
    }
    *mat = malloc(sizeof(struct matrix));
    if (*mat == NULL) {
      return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->data = calloc(rows * cols , sizeof(double));
    if ((*mat)->data == NULL) {
      return -1;
    }
    (*mat)->ref_cnt = 1;
    (*mat)->parent = NULL;
    return 0;
}

void print_matrix(matrix *mat) {
  printf("%d ", mat->rows);
  printf("%d\n", mat->cols);
  for (int i = 0; i < mat->rows; i++) {
    for (int j = 0; j < mat->cols; j++) {
      printf("%f ", mat->data[i * mat->cols + j]);
    }
    printf("%c", '\n');
  }
  printf("%c", '\n');
  printf("%c", '\n');
}


/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`.
 * You should return -1 if either `rows` or `cols` or both are non-positive or if any
 * call to allocate memory in this function fails. Return 0 upon success.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows <= 0 || cols <= 0) {
      return -1;
    }
    *mat = malloc(sizeof(struct matrix));
    if (*mat == NULL) {
      return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->data = from->data + offset;
    (*mat)->parent = from;
    (*mat)->parent->ref_cnt++;
    return 0;
}

/*
 * This function frees the matrix struct pointed to by `mat`. However, you need to make sure that
 * you only free the data if `mat` is not a slice and has no existing slices, or if `mat` is the
 * last existing slice of its parent matrix and its parent matrix has no other references.
 * You cannot assume that mat is not NULL.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat != NULL) {
        if (mat->parent == NULL) {
          if (mat->ref_cnt == 1) {
            free(mat->data);
            free(mat);
            return;
          } else if (mat->ref_cnt > 1) {
            mat->ref_cnt--;
            return;
          }
        } else {
          if (mat->parent->ref_cnt == 1) {
            free(mat->parent->data);
            free(mat->parent);
            free(mat);
            return;
          } else if (mat->parent->ref_cnt > 1) {
            mat->parent->ref_cnt--;
            free(mat);
            return;
          }
        }
    }
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    return mat->data[(row * mat->cols) + col];
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    mat->data[(row * mat->cols) + col] = val;
}

/*
 * Sets all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    for (int i = 0; i < (mat->rows * mat->cols); i++) {
      mat->data[i] = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if ((mat1->rows != mat2->rows) || (mat1->cols != mat2->cols)) {
      return 2;
    }
    if ((result->rows != mat1->rows) || (result->cols != mat1->cols)) {
      return 2;
    }
    if (mat1->cols <= 0 || mat2->cols <= 0 || mat1->rows <= 0 || mat2->rows <= 0) {
      return 1;
    }
    if ((mat1->rows * mat2->cols) >= 16) {
      #pragma omp parallel for
      for (int i = 0; i < ((mat1->rows * mat2->cols) / 16) * 16; i+=16) {

        __m256d _MAT1VALS1 = _mm256_loadu_pd((mat1->data + i));
        __m256d _MAT1VALS2 = _mm256_loadu_pd((mat1->data + i + 4));
        __m256d _MAT1VALS3 = _mm256_loadu_pd((mat1->data + i + 8));
        __m256d _MAT1VALS4 = _mm256_loadu_pd((mat1->data + i + 12));

        __m256d _MAT2VALS1 = _mm256_loadu_pd((mat2->data + i));
        __m256d _MAT2VALS2 = _mm256_loadu_pd((mat2->data + i + 4));
        __m256d _MAT2VALS3 = _mm256_loadu_pd((mat2->data + i + 8));
        __m256d _MAT2VALS4 = _mm256_loadu_pd((mat2->data + i + 12));

        __m256d _RESULTVALS1 = _mm256_add_pd(_MAT1VALS1, _MAT2VALS1);
        __m256d _RESULTVALS2 = _mm256_add_pd(_MAT1VALS2, _MAT2VALS2);
        __m256d _RESULTVALS3 = _mm256_add_pd(_MAT1VALS3, _MAT2VALS3);
        __m256d _RESULTVALS4 = _mm256_add_pd(_MAT1VALS4, _MAT2VALS4);

        _mm256_storeu_pd((result->data + i), _RESULTVALS1);
        _mm256_storeu_pd((result->data + i + 4), _RESULTVALS2);
        _mm256_storeu_pd((result->data + i + 8), _RESULTVALS3);
        _mm256_storeu_pd((result->data + i + 12), _RESULTVALS4);
      }
      for (int i = (mat1->rows * mat2->cols) / 16 * 16; i < (mat1->rows * mat2->cols); i++) {
        result->data[i] = mat1->data[i] + mat2->data[i];
      }
    } else {
      for (int i = 0; i <= (mat1->rows * mat2->cols); i++) {
        result->data[i] = mat1->data[i] + mat2->data[i];
      }
    }
    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if ((mat1->rows != mat2->rows) || (mat1->cols != mat2->cols)) {
      return 1;
    }
    if ((result->rows != mat1->rows) || (result->cols != mat1->cols)) {
      return 2;
    }
    if (mat1->cols <= 0 || mat2->cols <= 0 || mat1->rows <= 0 || mat2->rows <= 0) {
      return 1;
    }
    if ((mat1->rows * mat2->cols) >= 16) {
      #pragma omp parallel for
      for (int i = 0; i < ((mat1->rows * mat2->cols) / 16) * 16; i+=16) {

        __m256d _MAT1VALS1 = _mm256_loadu_pd((mat1->data + i));
        __m256d _MAT1VALS2 = _mm256_loadu_pd((mat1->data + i + 4));
        __m256d _MAT1VALS3 = _mm256_loadu_pd((mat1->data + i + 8));
        __m256d _MAT1VALS4 = _mm256_loadu_pd((mat1->data + i + 12));

        __m256d _MAT2VALS1 = _mm256_loadu_pd((mat2->data + i));
        __m256d _MAT2VALS2 = _mm256_loadu_pd((mat2->data + i + 4));
        __m256d _MAT2VALS3 = _mm256_loadu_pd((mat2->data + i + 8));
        __m256d _MAT2VALS4 = _mm256_loadu_pd((mat2->data + i + 12));

        __m256d _RESULTVALS1 = _mm256_sub_pd(_MAT1VALS1, _MAT2VALS1);
        __m256d _RESULTVALS2 = _mm256_sub_pd(_MAT1VALS2, _MAT2VALS2);
        __m256d _RESULTVALS3 = _mm256_sub_pd(_MAT1VALS3, _MAT2VALS3);
        __m256d _RESULTVALS4 = _mm256_sub_pd(_MAT1VALS4, _MAT2VALS4);

        _mm256_storeu_pd((result->data + i), _RESULTVALS1);
        _mm256_storeu_pd((result->data + i + 4), _RESULTVALS2);
        _mm256_storeu_pd((result->data + i + 8), _RESULTVALS3);
        _mm256_storeu_pd((result->data + i + 12), _RESULTVALS4);
      }
      for (int i = (((mat1->rows * mat2->cols) / 16) * 16); i < (mat1->rows * mat2->cols); i++) {
        result->data[i] = mat1->data[i] - mat2->data[i];
      }
    } else {
      for (int i = 0; i <= (mat1->rows * mat2->cols); i++) {
        result->data[i] = mat1->data[i] - mat2->data[i];
      }
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if ((mat1->cols != mat2->rows)) {
      return 1;
    }
    if ((result->rows != mat1->rows) || (result->cols != mat2->cols)) {
      return 2;
    }
    if (mat1->cols <= 0 || mat2->cols <= 0 || mat1->rows <= 0 || mat2->rows <= 0) {
      return 1;
    }
    int i, j, k;
    if (mat2->cols >= 16 && result->rows >= 16) {
      double tail[4];
      double tail1[4];
      double tail2[4];
      double tail3[4];
      double tail4[4];
      double tail5[4];
      double tail6[4];
      double tail7[4];
      double tail8[4];
      double tail9[4];
      double tail10[4];
      double tail11[4];
      double tail12[4];
      double tail13[4];
      double tail14[4];
      double tail15[4];
      double tail16[4];
      double tail17[4];
      double tail18[4];
      double tail19[4];
      double tail20[4];
      double tail21[4];
      double tail22[4];
      double tail23[4];
      double tail24[4];
      double tail25[4];
      double tail26[4];
      double tail27[4];
      double tail28[4];
      double tail29[4];
      double tail30[4];
      double tail31[4];
      #pragma omp parallel for private (k, tail, tail1, tail2, tail3, tail4, tail5, tail6, tail7, tail8, tail9, tail10, tail11, tail12, tail13, tail14, tail15, tail16, tail17, tail18, tail19, tail20, tail21, tail22, tail23, tail24, tail25, tail26, tail27, tail28, tail29, tail30, tail31)
      for (i = 0; i < result->rows / 32 * 32; i+=32) {
        for (j = 0; j < result->cols; j++) {
          tail[0] = tail[1] = tail[2] = tail[3] = 0;
          tail1[0] = tail1[1] = tail1[2] = tail1[3] = 0;
          tail2[0] = tail2[1] = tail2[2] = tail2[3] = 0;
          tail3[0] = tail3[1] = tail3[2] = tail3[3] = 0;
          tail4[0] = tail4[1] = tail4[2] = tail4[3] = 0;
          tail5[0] = tail5[1] = tail5[2] = tail5[3] = 0;
          tail6[0] = tail6[1] = tail6[2] = tail6[3] = 0;
          tail7[0] = tail7[1] = tail7[2] = tail7[3] = 0;
          tail8[0] = tail8[1] = tail8[2] = tail8[3] = 0;
          tail9[0] = tail9[1] = tail9[2] = tail9[3] = 0;
          tail10[0] = tail10[1] = tail10[2] = tail10[3] = 0;
          tail11[0] = tail11[1] = tail11[2] = tail11[3] = 0;
          tail12[0] = tail12[1] = tail12[2] = tail12[3] = 0;
          tail13[0] = tail13[1] = tail13[2] = tail13[3] = 0;
          tail14[0] = tail14[1] = tail14[2] = tail14[3] = 0;
          tail15[0] = tail15[1] = tail15[2] = tail15[3] = 0;
          tail16[0] = tail16[1] = tail16[2] = tail16[3] = 0;
          tail17[0] = tail17[1] = tail17[2] = tail17[3] = 0;
          tail18[0] = tail18[1] = tail18[2] = tail18[3] = 0;
          tail19[0] = tail19[1] = tail19[2] = tail19[3] = 0;
          tail20[0] = tail20[1] = tail20[2] = tail20[3] = 0;
          tail21[0] = tail21[1] = tail21[2] = tail21[3] = 0;
          tail22[0] = tail22[1] = tail22[2] = tail22[3] = 0;
          tail23[0] = tail23[1] = tail23[2] = tail23[3] = 0;
          tail24[0] = tail24[1] = tail24[2] = tail24[3] = 0;
          tail25[0] = tail25[1] = tail25[2] = tail25[3] = 0;
          tail26[0] = tail26[1] = tail26[2] = tail26[3] = 0;
          tail27[0] = tail27[1] = tail27[2] = tail27[3] = 0;
          tail28[0] = tail28[1] = tail28[2] = tail28[3] = 0;
          tail29[0] = tail29[1] = tail29[2] = tail29[3] = 0;
          tail30[0] = tail30[1] = tail30[2] = tail30[3] = 0;
          tail31[0] = tail31[1] = tail31[2] = tail31[3] = 0;
          for (k = 0; k < mat1->cols / 4 * 4; k += 4) {
            _mm256_storeu_pd((tail),
            _mm256_add_pd(*(__m256d *) tail,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + (i * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
          ))));
            _mm256_storeu_pd((tail1),
            _mm256_add_pd(*(__m256d *) tail1,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 1) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail2),
            _mm256_add_pd(*(__m256d *) tail2,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 2) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail3),
            _mm256_add_pd(*(__m256d *) tail3,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 3) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail4),
            _mm256_add_pd(*(__m256d *) tail4,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 4) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail5),
            _mm256_add_pd(*(__m256d *) tail5,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 5) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail6),
            _mm256_add_pd(*(__m256d *) tail6,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 6) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail7),
            _mm256_add_pd(*(__m256d *) tail7,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 7) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail8),
            _mm256_add_pd(*(__m256d *) tail8,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 8) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail9),
            _mm256_add_pd(*(__m256d *) tail9,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 9) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail10),
            _mm256_add_pd(*(__m256d *) tail10,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 10) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail11),
            _mm256_add_pd(*(__m256d *) tail11,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 11) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail12),
            _mm256_add_pd(*(__m256d *) tail12,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 12) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail13),
            _mm256_add_pd(*(__m256d *) tail13,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 13) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail14),
            _mm256_add_pd(*(__m256d *) tail14,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 14) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail15),
            _mm256_add_pd(*(__m256d *) tail15,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 15) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));

            /*******/

            _mm256_storeu_pd((tail16),
            _mm256_add_pd(*(__m256d *) tail16,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 16) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
          ))));
            _mm256_storeu_pd((tail17),
            _mm256_add_pd(*(__m256d *) tail17,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 17) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail18),
            _mm256_add_pd(*(__m256d *) tail18,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 18) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail19),
            _mm256_add_pd(*(__m256d *) tail19,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 19) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail20),
            _mm256_add_pd(*(__m256d *) tail20,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 20) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail21),
            _mm256_add_pd(*(__m256d *) tail21,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 21) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail22),
            _mm256_add_pd(*(__m256d *) tail22,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 22) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail23),
            _mm256_add_pd(*(__m256d *) tail23,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 23) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail24),
            _mm256_add_pd(*(__m256d *) tail24,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 24) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail25),
            _mm256_add_pd(*(__m256d *) tail25,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 25) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail26),
            _mm256_add_pd(*(__m256d *) tail26,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 26) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail27),
            _mm256_add_pd(*(__m256d *) tail27,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 27) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail28),
            _mm256_add_pd(*(__m256d *) tail28,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 28) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail29),
            _mm256_add_pd(*(__m256d *) tail29,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 29) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail30),
            _mm256_add_pd(*(__m256d *) tail30,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 30) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
            _mm256_storeu_pd((tail31),
            _mm256_add_pd(*(__m256d *) tail31,
            _mm256_mul_pd(*(__m256d *) (mat1->data + k + ((i + 31) * mat1->cols)),
            _mm256_set_pd(
              mat2->data[j + ((k + 3) * mat2->cols)],
              mat2->data[j + ((k + 2) * mat2->cols)],
              mat2->data[j + ((k + 1) * mat2->cols)],
              mat2->data[j + (k * mat2->cols)]
            ))));
          }
          result->data[j + i * result->cols] += tail[0] + tail[1] + tail[2] + tail[3];
          result->data[j + (i + 1) * result->cols] += tail1[0] + tail1[1] + tail1[2] + tail1[3];
          result->data[j + (i + 2) * result->cols] += tail2[0] + tail2[1] + tail2[2] + tail2[3];
          result->data[j + (i + 3) * result->cols] += tail3[0] + tail3[1] + tail3[2] + tail3[3];
          result->data[j + (i + 4) * result->cols] += tail4[0] + tail4[1] + tail4[2] + tail4[3];
          result->data[j + (i + 5) * result->cols] += tail5[0] + tail5[1] + tail5[2] + tail5[3];
          result->data[j + (i + 6) * result->cols] += tail6[0] + tail6[1] + tail6[2] + tail6[3];
          result->data[j + (i + 7) * result->cols] += tail7[0] + tail7[1] + tail7[2] + tail7[3];
          result->data[j + (i + 8) * result->cols] += tail8[0] + tail8[1] + tail8[2] + tail8[3];
          result->data[j + (i + 9) * result->cols] += tail9[0] + tail9[1] + tail9[2] + tail9[3];
          result->data[j + (i + 10) * result->cols] += tail10[0] + tail10[1] + tail10[2] + tail10[3];
          result->data[j + (i + 11) * result->cols] += tail11[0] + tail11[1] + tail11[2] + tail11[3];
          result->data[j + (i + 12) * result->cols] += tail12[0] + tail12[1] + tail12[2] + tail12[3];
          result->data[j + (i + 13) * result->cols] += tail13[0] + tail13[1] + tail13[2] + tail13[3];
          result->data[j + (i + 14) * result->cols] += tail14[0] + tail14[1] + tail14[2] + tail14[3];
          result->data[j + (i + 15) * result->cols] += tail15[0] + tail15[1] + tail15[2] + tail15[3];


          result->data[j + (i + 16) * result->cols] += tail16[0] + tail16[1] + tail16[2] + tail16[3];
          result->data[j + (i + 17) * result->cols] += tail17[0] + tail17[1] + tail17[2] + tail17[3];
          result->data[j + (i + 18) * result->cols] += tail18[0] + tail18[1] + tail18[2] + tail18[3];
          result->data[j + (i + 19) * result->cols] += tail19[0] + tail19[1] + tail19[2] + tail19[3];
          result->data[j + (i + 20) * result->cols] += tail20[0] + tail20[1] + tail20[2] + tail20[3];
          result->data[j + (i + 21) * result->cols] += tail21[0] + tail21[1] + tail21[2] + tail21[3];
          result->data[j + (i + 22) * result->cols] += tail22[0] + tail22[1] + tail22[2] + tail22[3];
          result->data[j + (i + 23) * result->cols] += tail23[0] + tail23[1] + tail23[2] + tail23[3];
          result->data[j + (i + 24) * result->cols] += tail24[0] + tail24[1] + tail24[2] + tail24[3];
          result->data[j + (i + 25) * result->cols] += tail25[0] + tail25[1] + tail25[2] + tail25[3];
          result->data[j + (i + 26) * result->cols] += tail26[0] + tail26[1] + tail26[2] + tail26[3];
          result->data[j + (i + 27) * result->cols] += tail27[0] + tail27[1] + tail27[2] + tail27[3];
          result->data[j + (i + 28) * result->cols] += tail28[0] + tail28[1] + tail28[2] + tail28[3];
          result->data[j + (i + 29) * result->cols] += tail29[0] + tail29[1] + tail29[2] + tail29[3];
          result->data[j + (i + 30) * result->cols] += tail30[0] + tail30[1] + tail30[2] + tail30[3];
          result->data[j + (i + 31) * result->cols] += tail31[0] + tail31[1] + tail31[2] + tail31[3];
          for (k = mat1->cols / 4 * 4; k < mat1->cols; k++) {
            result->data[j + (i * result->cols)] += mat1->data[k + i * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 1) * result->cols)] += mat1->data[k + (i + 1) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 2) * result->cols)] += mat1->data[k + (i + 2) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 3) * result->cols)] += mat1->data[k + (i + 3) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 4) * result->cols)] += mat1->data[k + (i + 4) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 5) * result->cols)] += mat1->data[k + (i + 5) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 6) * result->cols)] += mat1->data[k + (i + 6) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 7) * result->cols)] += mat1->data[k + (i + 7) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 8) * result->cols)] += mat1->data[k + (i + 8) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 9) * result->cols)] += mat1->data[k + (i + 9) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 10) * result->cols)] += mat1->data[k + (i + 10) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 11) * result->cols)] += mat1->data[k + (i + 11) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 12) * result->cols)] += mat1->data[k + (i + 12) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 13) * result->cols)] += mat1->data[k + (i + 13) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 14) * result->cols)] += mat1->data[k + (i + 14) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 15) * result->cols)] += mat1->data[k + (i + 15) * mat1->cols] * mat2->data[j + k * mat2->cols];

            result->data[j + ((i + 16) * result->cols)] += mat1->data[k + (i + 16) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 17) * result->cols)] += mat1->data[k + (i + 17) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 18) * result->cols)] += mat1->data[k + (i + 18) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 19) * result->cols)] += mat1->data[k + (i + 19) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 20) * result->cols)] += mat1->data[k + (i + 20) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 21) * result->cols)] += mat1->data[k + (i + 21) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 22) * result->cols)] += mat1->data[k + (i + 22) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 23) * result->cols)] += mat1->data[k + (i + 23) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 24) * result->cols)] += mat1->data[k + (i + 24) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 25) * result->cols)] += mat1->data[k + (i + 25) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 26) * result->cols)] += mat1->data[k + (i + 26) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 27) * result->cols)] += mat1->data[k + (i + 27) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 28) * result->cols)] += mat1->data[k + (i + 28) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 29) * result->cols)] += mat1->data[k + (i + 29) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 30) * result->cols)] += mat1->data[k + (i + 30) * mat1->cols] * mat2->data[j + k * mat2->cols];
            result->data[j + ((i + 31) * result->cols)] += mat1->data[k + (i + 31) * mat1->cols] * mat2->data[j + k * mat2->cols];
          }
          }
      }
      for (i = result->rows / 32 * 32; i < result->rows; i++) {
        for (j = 0; j < result->cols; j++) {
          for (k = 0; k < mat1->cols; k++) {
            result->data[j + i * result->cols] += mat1->data[k + i * mat1->cols] * mat2->data[j + k * mat2->cols];
          }
        }
      }
    } else {
      for (i = 0; i < result->rows; i++)
        for (j = 0; j < result->cols; j++) {
          for (k = 0; k < mat1->cols; k++)
            result->data[j + i * result->cols] += mat1->data[k + i * mat1->cols] * mat2->data[j + k * mat2->cols];
          }
    }

    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    if ((pow <= 0) || (mat->rows != mat->cols)) {
        return 1;
    }
    if ((result->rows != mat->rows) || (result->cols != mat->cols)) {
      return 2;
    }
    if (mat->cols <= 0 || mat->rows <= 0) {
      return 1;
    }
    if (pow == 1) {
      for (int i = 0; i < (mat->rows * mat->cols) / 4 * 4; i+=4) {
        result->data[i] = mat->data[i];
        result->data[i + 1] = mat->data[i + 1];
        result->data[i + 2] = mat->data[i + 2];
        result->data[i + 3] = mat->data[i + 3];
      }
      for (int i = (mat->rows * mat->cols) / 4 * 4; i < mat->rows * mat->cols; i++) {
        result->data[i] = mat->data[i];
      }
      return 0;
    } else if (pow == 2) {
      mul_matrix(result, mat, mat);
      return 0;
    } else {
      matrix *res = NULL;
      allocate_matrix(&res, mat->rows, mat->cols);
      pow_matrix(res, mat, pow - 1);
      fill_matrix(result, 0);
      mul_matrix(result, res, mat);
      deallocate_matrix(res);
    }
    return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if ((result->rows != mat->rows) || (result->cols != mat->cols)) {
        return 2;
    }
    if (mat->cols <= 0 || mat->rows <= 0) {
      return 1;
    }
    if ((mat->rows * mat->cols >= 16)) {
      __m256d _ZEROES = _mm256_set1_pd(0);
      #pragma omp parallel for
      for (int i = 0; i < ((mat->rows * mat->cols) / 16) * 16; i+=16) {

        __m256d _MATVALS1 = _mm256_loadu_pd((mat->data + i));
        __m256d _MATVALS2 = _mm256_loadu_pd((mat->data + i + 4));
        __m256d _MATVALS3 = _mm256_loadu_pd((mat->data + i + 8));
        __m256d _MATVALS4 = _mm256_loadu_pd((mat->data + i + 12));

        __m256d _RESULTVALS1 = _mm256_sub_pd(_ZEROES, _MATVALS1);
        __m256d _RESULTVALS2 = _mm256_sub_pd(_ZEROES, _MATVALS2);
        __m256d _RESULTVALS3 = _mm256_sub_pd(_ZEROES, _MATVALS3);
        __m256d _RESULTVALS4 = _mm256_sub_pd(_ZEROES, _MATVALS4);

        _mm256_storeu_pd((result->data + i), _RESULTVALS1);
        _mm256_storeu_pd((result->data + i + 4), _RESULTVALS2);
        _mm256_storeu_pd((result->data + i + 8), _RESULTVALS3);
        _mm256_storeu_pd((result->data + i + 12), _RESULTVALS4);
      }
      for (int i = ((mat->rows * mat->cols) / 16) * 16; i < (mat->rows * mat->cols); i++) {
        result->data[i] = mat->data[i] * -1;
      }
    } else {
      for (int i = 0; i < (mat->rows * mat->cols); i++) {
        result->data[i] = mat->data[i] * -1;
      }
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if ((result->rows != mat->rows) || (result->cols != mat->cols)) {
        return 2;
    }
    if (mat->cols <= 0 || mat->rows <= 0) {
      return 1;
    }
    #pragma omp parallel for
    for (int i = 0; i < ((mat->rows * mat->cols) / 4) * 4; i+=4) {

      if (mat->data[i] < 0) {
        result->data[i] = mat->data[i] * -1;
      } else {
        result->data[i] = mat->data[i];
      }

      if (mat->data[i + 1] < 0) {
        result->data[i + 1] = mat->data[i + 1] * -1;
      } else {
        result->data[i + 1] = mat->data[i + 1];
      }

      if (mat->data[i + 2] < 0) {
        result->data[i + 2] = mat->data[i + 2] * -1;
      } else {
        result->data[i + 2] = mat->data[i + 2];
      }

      if (mat->data[i + 3] < 0) {
        result->data[i + 3] = mat->data[i + 3] * -1;
      } else {
        result->data[i + 3] = mat->data[i + 3];
      }

    }

    for (int i = ((mat->rows * mat->cols) / 4) * 4; i < (mat->rows * mat->cols); i++) {
      if (mat->data[i] < 0) {
        result->data[i] = mat->data[i] * -1;
      } else {
        result->data[i] = mat->data[i];
      }
    }
    return 0;
}
