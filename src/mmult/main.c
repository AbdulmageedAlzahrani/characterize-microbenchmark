/* main.c
 *
 * Author: Khalid Al-Hawaj
 * Date  : 12 Nov. 2023
 *
 * This file is structured to call different implementation of the same
 * algorithm/microbenchmark. The file will allocate 3 output arrays one
 * for: scalar naive impl, scalar opt impl, vectorized impl. As it stands
 * the file will allocate and initialize with random data one input array
 * of type 'byte'. To check correctness, the file allocate a 'ref' array;
 * to calculate this 'ref' array, the file will invoke a ref_impl, which
 * is supposed to be functionally correct and act as a reference for
 * the functionality. The file also adds a guard word at the end of the
 * output arrays to check for buffer overruns.
 *
 * The file will invoke each implementation n number of times. It will
 * record the runtime of _each_ invocation through the following Linux
 * API:
 *    clock_gettime(), with the clk_id set to CLOCK_MONOTONIC
 * Then, the file will calculate the standard deviation and calculate
 * an outlier-free average by excluding runtimes that are larger than
 * 2 standard deviation of the original average.
 */

/* Set features         */
#define _GNU_SOURCE

/* Standard C includes  */
/*  -> Standard Library */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>


#include <fstream>
#include <sstream>
/*  -> Scheduling       */
#include <sched.h>
/*  -> Types            */
#include <stdbool.h>
#include <inttypes.h>
/*  -> Runtimes         */
#include <time.h>
#include <unistd.h>
#include <errno.h>

/* Include all implementations declarations */
#include "impl/ref.h"
#include "impl/naive.h"
#include "impl/opt.h"
#include "impl/vec.h"
#include "impl/para.h"

/* Include common headers */
#include "common/types.h"
#include "common/macros.h"

/* Include application-specific headers */
#include "include/types_m.h"

const int ROWS1 =  20 ;
const int COLS1 =  20 ;
const int ROWS2 =  20 ;
const int COLS2 =  20 ;

void initialize_matrix(float* matrix,int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // matrix[i * cols + j] = ((float)rand() / (float)RAND_MAX) * 1000.0f; 
            matrix[i * cols + j] = rand()%10+1;

        }
    }
}
static void load_matrix_from_csv(const char *filename,
  float *dest,
  int rows, int cols)
{
  
FILE* fp = fopen(filename, "r");

// print fp value
printf("fp: %p\n", fp);
if (!fp) {
fprintf(stderr, "Error: cannot open file %s\n", filename);
return;
}

// We expect each row in the CSV file to have `cols` float values
// separated by commas, e.g.:
// 0.1234,0.5678,...,0.9999
for (int r = 0; r < rows; r++) {
for (int c = 0; c < cols; c++) {
// Read a float followed by a comma (except maybe the last in line)
if (fscanf(fp, "%f%*c", &dest[r * cols + c]) != 1) {
fprintf(stderr, "Error reading float from %s at row %d col %d\n",
filename, r, c);
fclose(fp);
return;
}
}
}

fclose(fp);
}
void fillMatrices(byte *src1, 
  byte *src2, 
  byte *ref,
  int rows1, int cols1,
  int rows2, int cols2,
  int rowsC, int colsC)
{
// Cast the byte* pointers to float* to store float data.
float *fA = (float*)src1;
float *fB = (float*)src2;
float *fC = (float*)ref;

// Load A, B, C from CSV files created by Julia
load_matrix_from_csv("src/mmult/test_data/A.csv", fA, rows1, cols1);
load_matrix_from_csv("src/mmult/test_data/B.csv", fB, rows2, cols2);
load_matrix_from_csv("src/mmult/test_data/C.csv", fC, rowsC, colsC);
}



int main(int argc, char** argv)
{
  /* Set the buffer for printf to NULL */
  setbuf(stdout, NULL);

  /* Arguments */
  int nthreads = 1;
  int cpu      = 0;

  int nruns    = 10000;
  int nstdevs  = 3;

  /* Data */
  //int data_size = SIZE_DATA;
  int rows1 = ROWS1;
  int cols1 = COLS1;
  int rows2 = ROWS2;
  int cols2 = COLS2;

  /* Parse arguments */
  /* Function pointers */
  void* (*impl_scalar_naive_ptr)(void* args) = impl_scalar_naive;
  void* (*impl_scalar_opt_ptr  )(void* args) = impl_scalar_opt;
  void* (*impl_vector_ptr      )(void* args) = impl_vector;
  void* (*impl_parallel_ptr    )(void* args) = impl_parallel;

  /* Chosen */
  void* (*impl)(void* args) = NULL;
  const char* impl_str      = NULL;
  bool run_test = false;


  bool help = false;
  for (int i = 1; i < argc; i++) {
    /* Implementations */
    if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--impl") == 0) {
      assert (++i < argc);
      if (strcmp(argv[i], "naive") == 0) {
        impl = impl_scalar_naive_ptr; impl_str = "scalar_naive";
      } else if (strcmp(argv[i], "opt"  ) == 0) {
        impl = impl_scalar_opt_ptr  ; impl_str = "scalar_opt"  ;
      } else if (strcmp(argv[i], "vec"  ) == 0) {
        impl = impl_vector_ptr      ; impl_str = "vectorized"  ;
      } else if (strcmp(argv[i], "para" ) == 0) {
        impl = impl_parallel_ptr    ; impl_str = "parallelized";
      } else {
        impl = NULL                 ; impl_str = "unknown"     ;
      }

      continue;
    }

    /* Input/output data size */
    if (strcmp(argv[i], "-r1") == 0 || strcmp(argv[i], "--row1") == 0) {
      assert (++i < argc);
      rows1 = atoi(argv[i]);

      continue;
    }


    if (strcmp(argv[i], "-l1") == 0 || strcmp(argv[i], "--cols1") == 0) {
      assert (++i < argc);
      cols1 = atoi(argv[i]);

      continue;
    }
    if(strcmp(argv[i], "-r2") == 0 || strcmp(argv[i], "--row2") == 0) {
      assert (++i < argc);
      rows2 = atoi(argv[i]);

      continue;
    }
    if(strcmp(argv[i], "-l2") == 0 || strcmp(argv[i], "--cols2") == 0) {
      assert (++i < argc);
      cols2 = atoi(argv[i]);

      continue;
    }

    /* Run parameterization */
    if (strcmp(argv[i], "--nruns") == 0) {
      assert (++i < argc);
      nruns = atoi(argv[i]);

      continue;
    }

    if (strcmp(argv[i], "--nstdevs") == 0) {
      assert (++i < argc);
      nstdevs = atoi(argv[i]);

      continue;
    }

    /* Parallelization */
    if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nthreads") == 0) {
      assert (++i < argc);
      nthreads = atoi(argv[i]);

      continue;
    }

    if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cpu") == 0) {
      assert (++i < argc);
      cpu = atoi(argv[i]);

      continue;
    }

    /* Help */
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      help = true;

      continue;
    }

    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--test") == 0) {
      run_test = true;

      continue;
    }
    
  }

  if (help || impl == NULL) {
    if (!help) {
      if (impl_str != NULL) {
        printf("\n");
        printf("ERROR: Unknown \"%s\" implementation.\n", impl_str);
      } else {
        printf("\n");
        printf("ERROR: No implementation was chosen.\n");
      }
    }
    printf("\n");
    printf("Usage:\n");
    printf("  %s {-i | --impl} impl_str [Options]\n", argv[0]);
    printf("  \n");
    printf("  Required:\n");
    printf("    -i | --impl      Available implementations = {naive, opt, vec, para}\n");
    printf("    \n");
    printf("  Options:\n");
    printf("    -h | --help      Print this message\n");
    printf("    -n | --nthreads  Set number of threads available (default = %d)\n", nthreads);
    printf("    -c | --cpu       Set the main CPU for the program (default = %d)\n", cpu);
    printf("    -r1 | --rows      Number of rows in first input matrix (default = %d)\n", rows1);
    printf("    -l1 | --cols      Number of columns in first input matrix (default = %d)\n", cols1);
    printf("    -r2 | --rows      Number of rows in second input matrix (default = %d)\n", rows2);
    printf("    -l2 | --cols      Number of columns in second input matrix (default = %d)\n", cols2);
    printf("         --nruns     Number of runs to the implementation (default = %d)\n", nruns);
    printf("         --stdevs    Number of standard deviation to exclude outliers (default = %d)\n", nstdevs);
    printf("     -t | --test"    "run the pre-calculated test matrices\n");
    printf("\n");

    exit(help? 0 : 1);
  }

  if(cols1 != rows2) {
    printf("ERROR: The number of columns in the first matrix must be equal to the number of rows in the second matrix\n");
    exit(1);
  }

  /* Set our priority the highest */
  int nice_level = -20;

  printf("Setting up schedulers and affinity:\n");
  printf("  * Setting the niceness level:\n");
  do {
    errno = 0;
    printf("      -> trying niceness level = %d\n", nice_level);
    int __attribute__((unused)) ret = nice(nice_level);
  } while (errno != 0 && nice_level++);

  printf("    + Process has niceness level = %d\n", nice_level);

  /* If we are on an apple operating system, skip the scheduling  *
   * routine; Darwin does not support sched_set* system calls ... *
   *                                                              *
   * hawajkm: and here I was--thinking that MacOS is POSIX ...    *
   *          Silly me!                                           */
#if !defined(__APPLE__)
  /* Set scheduling to reduce context switching */
  /*    -> Set scheduling scheme                */
  printf("  * Setting up FIFO scheduling scheme and high priority ... ");
  pid_t pid    = 0;
  int   policy = SCHED_FIFO;
  struct sched_param param;

  param.sched_priority = sched_get_priority_max(policy);
  int res = sched_setscheduler(pid, policy, &param);
  if (res != 0) {
    printf("Failed\n");
  } else {
    printf("Succeeded\n");
  }

  /*    -> Set affinity                         */
  printf("  * Setting up scheduling affinity ... ");
  cpu_set_t cpumask;

  CPU_ZERO(&cpumask);
  for (int i = 0; i < nthreads; i++) {
    CPU_SET((cpu + i) % nthreads, &cpumask);
  }

  res = sched_setaffinity(pid, sizeof(cpumask), &cpumask);

  if (res != 0) {
    printf("Failed\n");
  } else {
    printf("Succeeded\n");
  }
#endif
  printf("\n");


  /* Statistics */
  __DECLARE_STATS(nruns, nstdevs);


  /* Initialize Rand */
  srand(0xdeadbeef);

  /* Datasets */
  /* Allocation and initialization */

  byte* src1   = __ALLOC_DATA     (byte, sizeof(float)*rows1*cols1);
  initialize_matrix((float*)src1,rows1,cols1); 
  byte* src2   = __ALLOC_DATA     (byte, sizeof(float)*rows2*cols2);
  initialize_matrix((float*)src2,rows2,cols2); 
  byte* ref   = __ALLOC_DATA     (byte, sizeof(float)*(rows1*cols2 + 1));
  byte* dest  = __ALLOC_DATA     (byte, sizeof(float)*(rows1*cols2 + 1));
  
  
  if (run_test) {
    printf("Running test matrices\n");
    rows1 = 5;
    cols1 = 5;
    rows2 = 5;
    cols2 = 7;

    fillMatrices(
      src1, src2, ref, 
      rows1, cols1, rows2, cols2, rows1, cols2);
  }

  // print rows and cols
  printf("Rows1: %d\n", rows1);
  printf("Cols1: %d\n", cols1);
  printf("Rows2: %d\n", rows2);
  printf("Cols2: %d\n", cols2);

  /* Setting a guards, which is 0xdeadcafe.
     The guard should not change or be touched. */

  __SET_GUARD(ref , sizeof(float)*rows1*cols2);
  __SET_GUARD(dest, sizeof(float)*rows1*cols2);


  /* Generate ref data */
  /* Arguments for the functions */
  args_t args_ref;



  args_ref.rows1     = rows1;
  args_ref.cols1     = cols1;
  args_ref.rows2     = rows2;
  args_ref.cols2     = cols2;
  args_ref.input1    = src1;
  args_ref.input2    = src2;
  args_ref.output    = ref;

  args_ref.cpu      = cpu;
  args_ref.nthreads = nthreads;





  /* Running the reference function */
  if (!run_test)
  impl_ref(&args_ref);


  /* Execute the requested implementation */
  /* Arguments for the function */
  args_t args;



  args.rows1     = rows1;
  args.cols1     = cols1;
  args.rows2     = rows2;
  args.cols2     = cols2;
  args.input1    = src1;
  args.input2    = src2;
  args.output    = dest;

  args.cpu      = cpu;
  args.nthreads = nthreads;



  /* Start execution */
  printf("Running \"%s\" implementation:\n", impl_str);

  printf("  * Invoking the implementation %d times .... ", num_runs);
  for (int i = 0; i < num_runs; i++) {
    __SET_START_TIME();
    for (int j = 0; j < 16; j++) {

      (*impl)(&args);
    }
    __SET_END_TIME();
    runtimes[i] = __CALC_RUNTIME() / 16;
  }
  printf("Finished\n");

  /* Verfication */
  printf("  * Verifying results .... ");
  bool match = __CHECK_FLOAT_MATCH(ref, dest, rows1*cols2-1, 0.1);
  bool guard = __CHECK_GUARD(     dest, sizeof(float)*rows1*cols2);

    // print src
  
  // FILE * fp_dump;
  // fp_dump = fopen("result.txt", "w");

  // if (fp_dump != NULL) {
  //    fprintf(fp_dump,"Src1:\n");
  // for (int i = 0; i < rows1; i++) {
  //   for (int j = 0; j < cols1; j++) {
  //     fprintf(fp_dump,"%f ", ((float*)src1)[i * cols1 + j]);
  //   }
  //   fprintf(fp_dump,"\n");
  // }
  // fprintf(fp_dump,"Src2:\n");
  // for (int i = 0; i < rows2; i++) {
  //   for (int j = 0; j < cols2; j++) {
  //     fprintf(fp_dump,"%f ", ((float*)src2)[i * cols2 + j]);
  //   }
  //   fprintf(fp_dump,"\n");
  // }

  // //print ref and dest
  // fprintf(fp_dump,"Ref:\n");
  // for (int i = 0; i < rows1; i++) {
  //   for (int j = 0; j < cols2; j++) {
  //     fprintf(fp_dump,"%f ", ((float*)ref)[i * cols2 + j]);
  //   }
  //   fprintf(fp_dump,"\n");
  // }
  // fprintf(fp_dump,"Dest:\n");
  // for (int i = 0; i < rows1; i++) {
  //   for (int j = 0; j < cols2; j++) {
  //     fprintf(fp_dump,"%f ", ((float*)dest)[i * cols2 + j]);
  //   }
  //   fprintf(fp_dump,"\n");
  // }

  

  //   fclose(fp_dump);

  // } 



  
 

  if (match && guard) {
    printf("Success\n");
  } else if (!match && guard) {
    printf("Fail, but no buffer overruns\n");
  } else if (match && !guard) {
    printf("Success, but failed buffer overruns check\n");
  } else if(!match && !guard) {
    printf("Failed, and failed buffer overruns check\n");
  }

  /* Running analytics */
  uint64_t min     = -1;
  uint64_t max     =  0;

  uint64_t avg     =  0;
  uint64_t avg_n   =  0;

  uint64_t std     =  0;
  uint64_t std_n   =  0;

  int      n_msked =  0;
  int      n_stats =  0;

  for (int i = 0; i < num_runs; i++)
    runtimes_mask[i] = true;

  printf("  * Running statistics:\n");
  do {
    n_stats++;
    printf("    + Starting statistics run number #%d:\n", n_stats);
    avg_n =  0;
    avg   =  0;

    /*   -> Calculate min, max, and avg */
    for (int i = 0; i < num_runs; i++) {
      if (runtimes_mask[i]) {
        if (runtimes[i] < min) {
          min = runtimes[i];
        }
        if (runtimes[i] > max) {
          max = runtimes[i];
        }
        avg += runtimes[i];
        avg_n += 1;
      }
    }
    avg = avg / avg_n;

    /*   -> Calculate standard deviation */
    std   =  0;
    std_n =  0;

    for (int i = 0; i < num_runs; i++) {
      if (runtimes_mask[i]) {
        std   += ((runtimes[i] - avg) *
                  (runtimes[i] - avg));
        std_n += 1;
      }
    }
    std = sqrt(std / std_n);

    /*   -> Calculate outlier-free average (mean) */
    n_msked = 0;
    for (int i = 0; i < num_runs; i++) {
      if (runtimes_mask[i]) {
        if (runtimes[i] > avg) {
          if ((runtimes[i] - avg) > (nstd * std)) {
            runtimes_mask[i] = false;
            n_msked += 1;
          }
        } else {
          if ((avg - runtimes[i]) > (nstd * std)) {
            runtimes_mask[i] = false;
            n_msked += 1;
          }
        }
      }
    }

    printf("      - Standard deviation = %" PRIu64 "\n", std);
    printf("      - Average = %" PRIu64 "\n", avg);
    printf("      - Number of active elements = %" PRIu64 "\n", avg_n);
    printf("      - Number of masked-off = %d\n", n_msked);
  } while (n_msked > 0);
  /* Display information */
  printf("  * Runtimes (%s): ", __PRINT_MATCH(match));
  printf(" %" PRIu64 " ns\n"  , avg                 );

  /* Dump */
  printf("  * Dumping runtime informations:\n");
  FILE * fp;
  char filename[256];
  strcpy(filename, impl_str);
  strcat(filename, "_runtimes.csv");
  printf("    - Filename: %s\n", filename);
  printf("    - Opening file .... ");
  fp = fopen(filename, "w");

  if (fp != NULL) {
    printf("Succeeded\n");
    printf("    - Writing runtimes ... ");
    fprintf(fp, "impl,%s", impl_str);

    fprintf(fp, "\n");
    fprintf(fp, "num_of_runs,%d", num_runs);

    fprintf(fp, "\n");
    fprintf(fp, "runtimes");
    for (int i = 0; i < num_runs; i++) {
      fprintf(fp, ", ");
      fprintf(fp, "%" PRIu64 "", runtimes[i]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "avg,%" PRIu64 "", avg);
    printf("Finished\n");
    printf("    - Closing file handle .... ");
    fclose(fp);
    printf("Finished\n");
  } else {
    printf("Failed\n");
  }
  printf("\n");

  /* Manage memory */
  free(src1);
  free(src2);
  free(dest);
  free(ref);

  /* Finished with statistics */
  __DESTROY_STATS();

  /* Done */
  return 0;
}
