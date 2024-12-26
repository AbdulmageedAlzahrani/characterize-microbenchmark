/* types.h
 *
 * Author: Khalid Al-Hawaj
 * Date  : 13 Nov. 2023
 * 
 * This file contains all required types decalartions.
*/

#ifndef __INCLUDE_TYPES_H_
#define __INCLUDE_TYPES_H_

typedef struct {
    const float* src0; // A
    const float* src1; // B
    float* dest; // C
    int rowsA;   // total number of rows in A
    int colsA;   // total number of cols in A (also rows in B)
    int colsB;   // total number of cols in B (also in C)
    int startRow;
    int endRow;
} ThreadData;

typedef struct {
  byte*   input1;
  byte*   input2;
  byte*   output;

  // size_t size;

  //size_t size;
  size_t rows1;
  size_t cols1;
  size_t rows2;
  size_t cols2;
  int     cpu;
  int     nthreads;
} args_t;

#endif //__INCLUDE_TYPES_H_
