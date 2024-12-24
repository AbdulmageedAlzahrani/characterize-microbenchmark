/* para.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <stdlib.h>
#include <math.h>

/* Include common headers */
#include "common/macros.h"
#include "common/types.h"
#include "string.h"

#include <ctype.h>
#include <stdio.h>

/* Include application-specific headers */
#include "include/types.h"
#include <bits/pthreadtypes.h>

/* Alternative Implementation */
#define inv_sqrt_2xPI 0.39894228040143270286

float CNDF(float InputX)
{
  int sign;

  float OutputX;
  float xInput;
  float xNPrimeofX;
  float expValues;
  float xK2;
  float xK2_2, xK2_3;
  float xK2_4, xK2_5;
  float xLocal, xLocal_1;
  float xLocal_2, xLocal_3;

  // Check for negative value of InputX
  if (InputX < 0.0)
  {
    InputX = -InputX;
    sign = 1;
  }
  elsenthreads
    sign = 0;

  xInput = InputX;

  // Compute NPrimeX term common to both four & six decimal accuracy calcs
  expValues = exp(-0.5f * InputX * InputX);
  xNPrimeofX = expValues;
  xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

  xK2 = 0.2316419 * xInput;
  xK2 = 1.0 + xK2;
  xK2 = 1.0 / xK2;
  xK2_2 = xK2 * xK2;
  xK2_3 = xK2_2 * xK2;
  xK2_4 = xK2_3 * xK2;
  xK2_5 = xK2_4 * xK2;

  xLocal_1 = xK2 * 0.319381530;
  xLocal_2 = xK2_2 * (-0.356563782);
  xLocal_3 = xK2_3 * 1.781477937;
  xLocal_2 = xLocal_2 + xLocal_3;
  xLocal_3 = xK2_4 * (-1.821255978);
  xLocal_2 = xLocal_2 + xLocal_3;
  xLocal_3 = xK2_5 * 1.330274429;
  xLocal_2 = xLocal_2 + xLocal_3;

  xLocal_1 = xLocal_2 + xLocal_1;
  xLocal = xLocal_1 * xNPrimeofX;
  xLocal = 1.0 - xLocal;

  OutputX = xLocal;

  if (sign)
  {
    OutputX = 1.0 - OutputX;
  }

  return OutputX;
}

float blackScholes(single_args_t *args)
{
  float OptionPrice;

  // local private working variables for the calculation
  float xStockPrice;
  float xStrikePrice;
  float xRiskFreeRate;
  float xVolatility;
  float xTime;
  float xSqrtTime;
  int otype;

  float logValues;
  float xLogTerm;
  float xD1;
  float xD2;
  float xPowerTerm;
  float xDen;
  float d1;
  float d2;
  float FutureValueX;
  float NofXd1;
  float NofXd2;
  float NegNofXd1;
  float NegNofXd2;

  xStockPrice = args->sptPrice;
  xStrikePrice = args->strike;
  xRiskFreeRate = args->rate;
  xVolatility = args->volatility;
  otype = args->otype;

  xTime = args->otime;
  xSqrtTime = sqrt(xTime);

  logValues = log(xStockPrice / xStrikePrice);

  xLogTerm = logValues;

  xPowerTerm = xVolatility * xVolatility;
  xPowerTerm = xPowerTerm * 0.5;

  xD1 = xRiskFreeRate + xPowerTerm;
  xD1 = xD1 * xTime;
  xD1 = xD1 + xLogTerm;

  xDen = xVolatility * xSqrtTime;
  xD1 = xD1 / xDen;
  xD2 = xD1 - xDen;

  d1 = xD1;
  d2 = xD2;

  NofXd1 = CNDF(d1);
  NofXd2 = CNDF(d2);

  FutureValueX = xStrikePrice * (exp(-(xRiskFreeRate) * (xTime)));
  if (otype == 0)
  {
    OptionPrice = (xStockPrice * NofXd1) - (FutureValueX * NofXd2);
  }
  else
  {
    NegNofXd1 = (1.0 - NofXd1);
    NegNofXd2 = (1.0 - NofXd2);
    OptionPrice = (FutureValueX * NegNofXd2) - (xStockPrice * NegNofXd1);
  }

  return OptionPrice;
}

// prepare the arguments for the parallel implementation based on the input arguments and index
single_args_t *prepare_args(args_t *args, int index)
{
  single_args_t *inputs = (single_args_t *)malloc(sizeof(args_t));
  if (inputs == NULL)
  {
    printf("Error: Cannot allocate memory for the arguments\n");
    exit(-2);
  }
  inputs->sptPrice = args->sptPrice[index];
  inputs->strike = args->strike[index];
  inputs->rate = args->rate[index];
  inputs->volatility = args->volatility[index];
  inputs->otime =  args->otime[index];
  inputs->otype = (tolower(args->otype[index]) == 'p')? 1 : 0;
  inputs->output = &args->output[index];

  return inputs;
}

void *impl_parallel(void *args)
{

  args_t *inputs = (args_t *)args;

  int nthreads = inputs->nthreads;
  int num_stocks = inputs->num_stocks;

  pthread_t threads[num_stocks];

  for (int i = 0; i < num_stocks; i++)
  {
    single_args_t *thread_args = prepare_args(inputs, i);
    pthread_create(&threads[i], NULL, blackScholes, (void *)thread_args);
  }

  for (int i = 0; i < num_stocks; i++)
  {
    pthread_join(threads[i], NULL);
  }

  return NULL;
}
