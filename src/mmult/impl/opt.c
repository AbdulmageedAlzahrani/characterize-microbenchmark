/* opt.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <stdlib.h>

/* Include common headers */
#include "common/macros.h"
#include "common/types.h"

/* Include application-specific headers */
#include "include/types_m.h"

/* Alternative Implementation */
// #pragma GCC push_options
// #pragma GCC optimize ("O1")
void* impl_scalar_opt(void* args)
{
    args_t *parsed_args = (args_t *)args;

    float*   dest = (      float*)(parsed_args->output);
  const float*   src0 = (const float*)(parsed_args->input1);
  const float*   src1 = (const float*)(parsed_args->input2);


  int rows1 = parsed_args->rows1;
  int cols1 = parsed_args->cols1;

  int rows2 = parsed_args->rows2;
  int cols2 = parsed_args->cols2;

int blockSize = 256; 
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            dest[i * cols2 + j] = 0.0f;
        }
    }
    for (int i = 0; i < rows1; i += blockSize) {
        for (int j = 0; j < cols2; j += blockSize) {
            for (int k = 0; k < cols1; k += blockSize) {
                for (int ii = i; ii < i + blockSize && ii < rows1; ii++) {
                    for (int jj = j; jj < j + blockSize && jj < cols2; jj++) {
                        float sum = 0.0f;
                        for (int kk = k; kk < k + blockSize && kk < cols1; kk++) {
                            sum += src0[ii * cols1 + kk] * src1[kk * cols2 + jj];
                        }
                        dest[ii * cols2 + jj] += sum;
                    }
                }
            }
        }
    }
    return NULL;
}
// #pragma GCC pop_options
