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
  size_t num_stocks;

  float* sptPrice  ;
  float* strike    ;
  float* rate      ;
  float* volatility;
  float* otime     ;
  char * otype     ;
  float* output    ;

  int    cpu;
  int    nthreads;
} args_t;

//create a structure to hold the arguments, so that we can pass them to the threads
typedef struct {
  

  float* sptPrice  ;
  float* strike    ;
  float* rate      ;
  float* volatility;
  float* otime     ;
  char*  otype     ;
  float *output    ;
  int start_index  ;
  int end_index    ;

} single_args_t;

#endif //__INCLUDE_TYPES_H_
