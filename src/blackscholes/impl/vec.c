/* vec.c
 *
 * Author:
 * Date  :
 *
 *  Description
 */

/* Standard C includes */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/macros.h"
#include "common/types.h"
#include "common/vmath.h"
#include "include/types.h"
#include "scalar.h"

#define inv_sqrt_2xPI 0.39894228040143270286

static inline __m256 CNDF_avx(__m256 InputX) {
  // Create a mask for InputX < 0
  __m256 zero = _mm256_setzero_ps();
  __m256 signMask = _mm256_cmp_ps(InputX, zero, _CMP_LT_OQ);

  // Replace NaNs in signMask
  __m256 nanMask = _mm256_cmp_ps(InputX, InputX, _CMP_UNORD_Q);
  signMask = _mm256_or_ps(signMask, nanMask);

  // Compute absolute value of InputX
  __m256 absInputX =
      _mm256_blendv_ps(InputX, _mm256_sub_ps(zero, InputX), signMask);

  InputX = absInputX;

  __m256 expValues = _mm256_mul_ps(InputX, InputX);
  expValues = _mm256_mul_ps(expValues, _mm256_set1_ps(-0.5f));
  expValues = _mm256_exp_ps(expValues);
  __m256 xNPrimeofX = _mm256_mul_ps(expValues, _mm256_set1_ps(inv_sqrt_2xPI));

  __m256 xK2 =
      _mm256_fmadd_ps(InputX, _mm256_set1_ps(0.2316419f), _mm256_set1_ps(1.0f));
  //  print the xK2 values
  xK2 = _mm256_div_ps(_mm256_set1_ps(1.0f), xK2);

  __m256 xK2_2 = _mm256_mul_ps(xK2, xK2);
  __m256 xK2_3 = _mm256_mul_ps(xK2_2, xK2);
  __m256 xK2_4 = _mm256_mul_ps(xK2_3, xK2);
  __m256 xK2_5 = _mm256_mul_ps(xK2_4, xK2);

  __m256 xLocal_1 = _mm256_mul_ps(xK2, _mm256_set1_ps(0.319381530f));
  __m256 xLocal_2 = _mm256_mul_ps(xK2_2, _mm256_set1_ps(-0.356563782f));
  __m256 xLocal_3 = _mm256_mul_ps(xK2_3, _mm256_set1_ps(1.781477937f));
  xLocal_2 = _mm256_add_ps(xLocal_2, xLocal_3);
  xLocal_3 = _mm256_mul_ps(xK2_4, _mm256_set1_ps(-1.821255978f));
  xLocal_2 = _mm256_add_ps(xLocal_2, xLocal_3);
  xLocal_3 = _mm256_mul_ps(xK2_5, _mm256_set1_ps(1.330274429f));
  xLocal_2 = _mm256_add_ps(xLocal_2, xLocal_3);

  xLocal_1 = _mm256_add_ps(xLocal_1, xLocal_2);

  __m256 xLocal = _mm256_mul_ps(xLocal_1, xNPrimeofX);
  xLocal = _mm256_sub_ps(_mm256_set1_ps(1.0f), xLocal);

  __m256 OutputX = _mm256_blendv_ps(
      xLocal, _mm256_sub_ps(_mm256_set1_ps(1.0f), xLocal), signMask);

  return OutputX;
}

void blackScholes_avx(float *sptprice, float *strike, float *rate,
                      float *volatility, float *otime, char *otype,
                      float *output, int num_stocks) {
  int i;
  for (i = 0; i <= num_stocks - 8; i += 8) {
    __m256 xStockPrice = _mm256_loadu_ps(&sptprice[i]);
    __m256 xStrikePrice = _mm256_loadu_ps(&strike[i]);
    __m256 xRiskFreeRate = _mm256_loadu_ps(&rate[i]);
    __m256 xVolatility = _mm256_loadu_ps(&volatility[i]);
    __m256 xTime = _mm256_loadu_ps(&otime[i]);

    __m256 xSqrtTime = _mm256_sqrt_ps(xTime);
    __m256 xLogTerm = _mm256_log_ps(_mm256_div_ps(xStockPrice, xStrikePrice));

    __m256 xPowerTerm = _mm256_mul_ps(xVolatility, xVolatility);
    xPowerTerm = _mm256_mul_ps(xPowerTerm, _mm256_set1_ps(0.5f));

    __m256 xD1 = _mm256_add_ps(xRiskFreeRate, xPowerTerm);
    xD1 = _mm256_mul_ps(xD1, xTime);
    xD1 = _mm256_add_ps(xD1, xLogTerm);

    __m256 xDen = _mm256_mul_ps(xVolatility, xSqrtTime);
    xD1 = _mm256_div_ps(xD1, xDen);
    __m256 xD2 = _mm256_sub_ps(xD1, xDen);
    float debugInput[8];

    __m256 NofXd1 = CNDF_avx(xD1);
    __m256 NofXd2 = CNDF_avx(xD2);

    __m256 FutureValueX = _mm256_mul_ps(
        xStrikePrice,
        _mm256_exp_ps(_mm256_mul_ps(_mm256_set1_ps(-1.0f),
                                    _mm256_mul_ps(xRiskFreeRate, xTime))));

    __m256i otypeInt =
        _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *)&otype[i]));
    __m256i pMask = _mm256_cmpeq_epi32(otypeInt, _mm256_set1_epi32((int)'P'));
    __m256 otypeMask = _mm256_castsi256_ps(pMask);

    __m256 CallPrice = _mm256_sub_ps(_mm256_mul_ps(xStockPrice, NofXd1),
                                     _mm256_mul_ps(FutureValueX, NofXd2));
    __m256 PutPrice = _mm256_sub_ps(
        _mm256_mul_ps(FutureValueX,
                      _mm256_sub_ps(_mm256_set1_ps(1.0f), NofXd2)),
        _mm256_mul_ps(xStockPrice,
                      _mm256_sub_ps(_mm256_set1_ps(1.0f), NofXd1)));

    __m256 OptionPrice = _mm256_blendv_ps(CallPrice, PutPrice, otypeMask);
    _mm256_storeu_ps(&output[i], OptionPrice);
  }

  // Handle remaining stocks
  for (; i < num_stocks; i++) {
    output[i] = blackScholes_s(sptprice[i], strike[i], rate[i], volatility[i],
                               otime[i], otype[i]);
  }
}

void *impl_vector(void *args) {
  args_t *inputs = (args_t *)args;
  blackScholes_avx(inputs->sptPrice, inputs->strike, inputs->rate,
                   inputs->volatility, inputs->otime, inputs->otype,
                   inputs->output, inputs->num_stocks);
  return NULL;
}
