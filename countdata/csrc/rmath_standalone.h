/*
 * Math backend for countdata Python wrapper.
 *
 * Two modes:
 *   1. USE_RMATH (recommended): Link against R's standalone libRmath.
 *      Gives identical results to R's countdata package.
 *      Compile: gcc -DUSE_RMATH -I$(R_INCLUDE) ... -lRmath -lm
 *
 *   2. Fallback: Standalone C99 implementations.
 *      Uses system lgamma + hand-coded digamma/trigamma/pnorm/pchisq.
 *      bb.test matches R; ibb.test may have small precision differences
 *      due to C99 lgamma vs R's lgammafn differing by 1-2 ULP.
 */
#ifndef RMATH_STANDALONE_H
#define RMATH_STANDALONE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 *  Mode 1: Use R's standalone math library (libRmath)
 * ================================================================ */
#ifdef USE_RMATH

#define MATHLIB_STANDALONE
#include <Rmath.h>
/* Rmath.h provides: lgammafn, digamma, trigamma, pnorm, pchisq */

/* ================================================================
 *  Mode 2: Standalone C99 implementations
 * ================================================================ */
#else

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static inline double lgammafn(double x) {
    return lgamma(x);
}

/*
 * digamma: recurrence to x >= 12, then 8-term Bernoulli asymptotic.
 * B_{2k}/(2k) coefficients for k=1..8.
 */
static double digamma(double x) {
    double result = 0.0;
    if (x != x) return x;
    if (x < 0.0) return digamma(1.0 - x) - M_PI / tan(M_PI * x);
    while (x < 12.0) { result -= 1.0 / x; x += 1.0; }
    double r = 1.0 / x, r2 = r * r, t = r2;
    result += log(x) - 0.5 * r;
    result -= (1.0/12.0) * t;       t *= r2;
    result += (1.0/120.0) * t;      t *= r2;
    result -= (1.0/252.0) * t;      t *= r2;
    result += (1.0/240.0) * t;      t *= r2;
    result -= (5.0/660.0) * t;      t *= r2;
    result += (691.0/32760.0) * t;   t *= r2;
    result -= (1.0/12.0) * t;       t *= r2;
    result += (3617.0/8160.0) * t;
    return result;
}

/*
 * trigamma: recurrence to x >= 12, then 8-term Bernoulli asymptotic.
 * B_{2k} coefficients for k=1..8.
 */
static double trigamma(double x) {
    double result = 0.0;
    if (x != x) return x;
    if (x < 0.0) {
        double s = sin(M_PI * x);
        return -(M_PI * M_PI) / (s * s) + trigamma(1.0 - x);
    }
    while (x < 12.0) { result += 1.0 / (x * x); x += 1.0; }
    double r = 1.0 / x, r2 = r * r, t = r2 * r;
    result += r + 0.5 * r2;
    result += (1.0/6.0) * t;        t *= r2;
    result -= (1.0/30.0) * t;       t *= r2;
    result += (1.0/42.0) * t;       t *= r2;
    result -= (1.0/30.0) * t;       t *= r2;
    result += (5.0/66.0) * t;       t *= r2;
    result -= (691.0/2730.0) * t;    t *= r2;
    result += (7.0/6.0) * t;        t *= r2;
    result -= (3617.0/510.0) * t;
    return result;
}

static double pnorm(double x, double mu, double sigma, int lower_tail, int log_p) {
    double q = (x - mu) / sigma;
    double p;
    if (fabs(q) > 38.0) p = (q > 0) ? 1.0 : 0.0;
    else p = 0.5 * erfc(-q * M_SQRT1_2);
    if (!lower_tail) p = 1.0 - p;
    if (log_p) p = log(p);
    return p;
}

static double pgamma_lower(double a, double x) {
    if (x <= 0.0) return 0.0;
    if (x < a + 1.0) {
        double term = 1.0 / a, sum = term, ap = a;
        for (int n = 0; n < 2000; n++) {
            ap += 1.0; term *= x / ap; sum += term;
            if (fabs(term) < fabs(sum) * 1e-16) break;
        }
        return sum * exp(-x + a * log(x) - lgamma(a));
    } else {
        double f = 1.0, c = 1.0, d = 1.0 / (x + 1.0 - a);
        f = d;
        for (int n = 1; n <= 2000; n++) {
            double an = -((double)n) * ((double)n - a);
            double bn = x + 2.0 * n + 1.0 - a;
            d = bn + an * d; if (fabs(d) < 1e-30) d = 1e-30; d = 1.0 / d;
            c = bn + an / c; if (fabs(c) < 1e-30) c = 1e-30;
            double delta = c * d; f *= delta;
            if (fabs(delta - 1.0) < 1e-16) break;
        }
        return 1.0 - exp(-x + a * log(x) - lgamma(a)) * f;
    }
}

static double pchisq(double x, double df, int lower_tail, int log_p) {
    double p;
    if (x <= 0.0) p = 0.0;
    else p = pgamma_lower(df / 2.0, x / 2.0);
    if (!lower_tail) p = 1.0 - p;
    if (log_p) p = log(p);
    return p;
}

#endif /* USE_RMATH */

/* ---- Common R compatibility stubs ---- */
#ifndef USE_RMATH
#define Rprintf printf
#define R_FlushConsole() fflush(stdout)
typedef int Rboolean;
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifdef __GNUC__
__attribute__((unused))
#endif
static Rboolean R_ToplevelExec(void (*fun)(void*), void* data) {
    (void)fun; (void)data; return TRUE;
}
static void R_CheckUserInterrupt(void) {}
#define LibExtern
#else
/* When using Rmath, these are already defined */
#define R_FlushConsole() fflush(stdout)
#define LibExtern
#endif

#endif /* RMATH_STANDALONE_H */