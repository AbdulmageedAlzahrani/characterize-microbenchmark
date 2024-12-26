/* vec.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>

/* Include common headers */
#include "common/macros.h"
#include "common/types.h"

/* Include application-specific headers */
#include "include/types_m.h"

/* Alternative Implementation */

void* impl_vector(void* args)
{
  args_t *parsed_args = (args_t *)args;

  float*   dest = (      float*)(parsed_args->output);
  const float*   src0 = (const float*)(parsed_args->input1);
  const float*   src1 = (const float*)(parsed_args->input2);
  // register       size_t size =              parsed_args->size / 4;

  int rows1 = parsed_args->rows1;
  int cols1 = parsed_args->cols1;

  int rows2 = parsed_args->rows2;
  int cols2 = parsed_args->cols2;

  const int max_vlen = 8;  // AVX2 processes 8 floats at a time
    // Perform matrix multiplication
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            __m256 sum = _mm256_setzero_ps();  // Initialize sum to 0
            dest[i * cols2 + j] = 0;
            int k = 0;

            // Perform dot product in chunks of 8
            for (; k <= cols1 - max_vlen; k += max_vlen) {
                // Load 8 elements from row A[i]
                __m256 a = _mm256_loadu_ps(&src0[i * cols1 + k]);

                // Load 8 elements from column B[k][j] (k-th row of B)
                __m256 b = _mm256_setr_ps(
                src1[(k + 0) * cols2 + j],
                src1[(k + 1) * cols2 + j],
                src1[(k + 2) * cols2 + j],
                src1[(k + 3) * cols2 + j],
                src1[(k + 4) * cols2 + j],
                src1[(k + 5) * cols2 + j],
                src1[(k + 6) * cols2 + j],
                src1[(k + 7) * cols2 + j]
            );     // Multiply and accumulate
                sum = _mm256_fmadd_ps(a, b, sum);  // sum += a * b
            }            // Sum up the elements of the vector 'sum' into a scalar
                __m128 lo = _mm256_castps256_ps128(sum);  // Lower half
                __m128 hi = _mm256_extractf128_ps(sum, 1); // Upper half
                lo = _mm_add_ps(lo, hi);  // Add halves
                lo = _mm_hadd_ps(lo, lo); // Horizontal add
                lo = _mm_hadd_ps(lo, lo);
                float result = _mm_cvtss_f32(lo);
            // Handle remaining elements (tail case)
            for (; k < cols1; k++) {
                result += src0[i * cols1 + k] * src1[k * cols2 + j];
            }
            // Store result in C[i][j]
            dest[i * cols2 + j] = result;
        }
    }

  return NULL;
}
// #pragma GCC pop_options
