/* para.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <stdlib.h>
#include <immintrin.h>
#include <pthread.h>

/* Include common headers */
#include "common/macros.h"
#include "common/types.h"

/* If we are on Darwin, include the compatibility header */
#if defined(__APPLE__)
#include "common/mach_pthread_compatibility.h"
#endif

/* Include application-specific headers */
#include "include/types_m.h"


/* Alternative Implementation */

void* thread_matmul(void* arg) 
{
    ThreadData* data = (ThreadData*)arg;
    
    const float* src0 = (const float*)data->src0;  // A
    const float* src1 = (const float*)data->src1;  // B
    float* dest = data->dest;  // C
    
    int rowsA = data->rowsA;
    int colsA = data->colsA;
    int colsB = data->colsB;
    
    int startRow = data->startRow;
    int endRow   = data->endRow;

    const int max_vlen = 8;  // AVX processes 8 floats at a time

    // Perform the multiplication only for [startRow .. endRow)
    for (int i = startRow; i < endRow; i++) {
        for (int j = 0; j < colsB; j++) {

            // Initialize an AVX accumulator
            __m256 sum = _mm256_setzero_ps();
            float result = 0.0f;
            int k = 0;

            // Perform dot product in chunks of 8
            for (; k <= colsA - max_vlen; k += max_vlen) {
                // Load 8 elements from row A[i] 
                __m256 a = _mm256_loadu_ps(&src0[i * colsA + k]);

                // Load 8 elements from column B[k][j] 
                // using ascending order with _mm256_setr_ps
                __m256 b = _mm256_setr_ps(
                    src1[(k + 0) * colsB + j],
                    src1[(k + 1) * colsB + j],
                    src1[(k + 2) * colsB + j],
                    src1[(k + 3) * colsB + j],
                    src1[(k + 4) * colsB + j],
                    src1[(k + 5) * colsB + j],
                    src1[(k + 6) * colsB + j],
                    src1[(k + 7) * colsB + j]
                );
                
                // Multiply and accumulate
                sum = _mm256_fmadd_ps(a, b, sum);
            }

            // Reduce the 8-lane AVX sum to a single float
            __m128 lo = _mm256_castps256_ps128(sum);
            __m128 hi = _mm256_extractf128_ps(sum, 1);
            lo = _mm_add_ps(lo, hi);
            lo = _mm_hadd_ps(lo, lo);
            lo = _mm_hadd_ps(lo, lo);
            result = _mm_cvtss_f32(lo);

            // Handle leftover elements (the tail case)
            for (; k < colsA; k++) {
                result += src0[i * colsA + k] * src1[k * colsB + j];
            }

            // Store in C
            dest[i * colsB + j] = result;
        }
    }

    // End of the thread function
    pthread_exit(NULL);
}

void* impl_parallel(void* args)
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

  int numThreads = 4;
  pthread_t* threads = (pthread_t*)malloc(numThreads * sizeof(pthread_t));
  ThreadData* tdata = (ThreadData*)malloc(numThreads * sizeof(ThreadData));

   // Divide the total rows among the threads
  int chunkSize = rows1 / numThreads;
  int remainder = rows1 % numThreads;

  int currentRow = 0;
    for (int t = 0; t < numThreads; t++) {
        int startRow = currentRow;
        // last chunk can take extra if there's a remainder
        int thisChunk = chunkSize + ((t < remainder) ? 1 : 0);
        int endRow = startRow + thisChunk;
        currentRow = endRow;

        // Fill out the struct
        tdata[t].src0    = src0;
        tdata[t].src1    = src1;
        tdata[t].dest    = dest;
        tdata[t].rowsA   = rows1;
        tdata[t].colsA   = cols1;
        tdata[t].colsB   = cols2;
        tdata[t].startRow = startRow;
        tdata[t].endRow   = endRow;

        // Launch
        pthread_create(&threads[t], NULL, thread_matmul, &tdata[t]);
    }

    // Join all threads
    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }

    // At this point, 'dest' should contain the product of A x B
    
    // Cleanup
    free(threads);
    free(tdata);


  
  return NULL;
}
