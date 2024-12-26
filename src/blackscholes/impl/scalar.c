/* scalar.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/* Include common headers */
#include "common/macros.h"
#include "common/types.h"

/* Include application-specific headers */
#include "include/types.h"
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
  else
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

float blackScholes(float sptprice, float strike, float rate, float volatility,
                   float otime, char otype)
{
  float OptionPrice;

  // local private working variables for the calculation
  float xStockPrice;
  float xStrikePrice;
  float xRiskFreeRate;
  float xVolatility;
  float xTime;
  float xSqrtTime;

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

  xStockPrice = sptprice;
  xStrikePrice = strike;
  xRiskFreeRate = rate;
  xVolatility = volatility;

  xTime = otime;
  xSqrtTime = sqrt(xTime);

  logValues = log(sptprice / strike);

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

  FutureValueX = strike * (exp(-(rate) * (otime)));
  otype = (tolower(otype) == 'p')? 1 : 0;

  if (otype == 0)
  {
    OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
  }
  else
  {
    NegNofXd1 = (1.0 - NofXd1);
    NegNofXd2 = (1.0 - NofXd2);
    OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
  }

  return OptionPrice;
}

/* Naive Implementation */
void *impl_scalar(void *args)
{
  args_t *inputs = (args_t *)args;
  int num_stocks = inputs->num_stocks;

  for (int i = 0; i < num_stocks; i++)
  {
    inputs->output[i] = blackScholes(inputs->sptPrice[i], inputs->strike[i], inputs->rate[i],
                                     inputs->volatility[i], inputs->otime[i], inputs->otype[i]);
  }

  return NULL;
}
