/* ref.c
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

/* Reference Implementation */
void* impl_ref(void* args)
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

  for(int i = 0; i < rows1; i++) {
    for(int j = 0; j < cols2; j++) {
      float temp = 0.0f;
      for(int k = 0; k < cols2; k++) {
        temp += src0[i * cols2 + k] * src1[k * cols2 + j];
      }
      dest[i * cols2 + j] = temp;
    }
  }
  return NULL;
}
