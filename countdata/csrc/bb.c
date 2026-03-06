/*
 * Beta-binomial test - standalone C for Python ctypes.
 * Original: Thang V. Pham, t.pham@amsterdamumc.nl
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "rmath_standalone.h"

#define SMAX (TYPE)1e6
#define SMIN (TYPE)1e-6
#define MAX_INNER 200
#define ALPHA  (TYPE)0.1
#define BETA (TYPE)0.7
#define SMALL (TYPE)1e-12
#define forint(i, a, b) for (int i=(int)(a); i<(int)(b); i++)
#define TYPE double

typedef struct _bb_t {
    TYPE* a; TYPE* ta;
    int N, M;
    int* g_size; int* g_ind;
    int theta_equal;
    TYPE _m1, _size;
    TYPE* _a; TYPE* _ta; TYPE* _m1_array;
    int comp;
    TYPE f, f0;
} bb_t;

TYPE opt_find_eta(TYPE* v, TYPE g, TYPE h, TYPE p, TYPE lb, TYPE ub) {
    if (fabs(h) < SMALL) h = (h > (TYPE)0 ? SMALL : -SMALL);
    *v = -g / h;
    if (g * (*v) > 0) *v = -(*v);
    TYPE eta = (TYPE)1.0;
    if (p + *v >= ub) eta = (ub - SMALL - p) / *v;
    if (p + *v <= lb) eta = (lb + SMALL - p) / *v;
    if (eta < (TYPE)0.0) eta = (TYPE)0.0;
    return eta;
}

TYPE minimize1d(TYPE(*fval)(TYPE*, TYPE*, bb_t*, TYPE), bb_t* x, TYPE p0, TYPE lb, TYPE ub) {
    TYPE p = p0, old_p = p;
    int cc = 0;
    forint(i, 0, MAX_INNER) {
        TYPE f, g, h;
        f = fval(&g, &h, x, p);
        TYPE v;
        TYPE eta = opt_find_eta(&v, g, h, p, lb, ub);
        TYPE lambda = g * v;
        TYPE fnew = fval(0, 0, x, p + eta * v);
        while ((fnew - f) > ALPHA * lambda * eta) {
            eta *= BETA;
            fnew = fval(0, 0, x, p + eta * v);
            if (eta < (TYPE)1e-15) { eta = (TYPE)0.0; fnew = f; break; }
        }
        p += eta * v;
        if (fabs(fnew - f) < (TYPE)1e-8 && fabs(p - old_p) < (TYPE)1e-8) cc++; else cc = 0;
        if (cc > 1) break;
        old_p = p;
    }
    return p;
}

void bb_simple_estimate_tm(bb_t* x, TYPE* alp, TYPE* bet) {
    TYPE p = 0.0, m2 = 0.0;
    forint(i, 0, x->_size) {
        TYPE tmp = x->_a[i] / x->_ta[i];
        p += tmp; m2 += tmp * tmp;
    }
    p /= (TYPE)(x->_size); m2 /= (TYPE)(x->_size);
    if (p < SMALL) { *alp = 1.0; *bet = 1e4; }
    else {
        TYPE s;
        if ((m2 - p * p) < SMALL) s = 1e4;
        else { s = (p - m2) / (m2 - p * p); if (s < SMIN) s = SMIN; if (s > SMAX) s = SMAX; }
        *alp = p * s; *bet = (1.0 - p) * s;
    }
}

TYPE fval_s(TYPE* g, TYPE* h, bb_t* x, TYPE s) {
    TYPE m2 = 1.0 - x->_m1;
    if (g) {
        TYPE f = 0.0; *g = 0.0; *h = 0.0;
        TYPE lg_a0 = lgammafn(s), di_a0 = digamma(s), tri_a0 = trigamma(s);
        TYPE a1 = s * x->_m1, a2 = s * m2;
        TYPE lg_a1 = lgammafn(a1), di_a1 = x->_m1 * digamma(a1), tri_a1 = x->_m1 * x->_m1 * trigamma(a1);
        TYPE lg_a2 = lgammafn(a2), di_a2 = m2 * digamma(a2), tri_a2 = m2 * m2 * trigamma(a2);
        forint(i, 0, x->_size) {
            TYPE term2 = x->_ta[i] + s;
            TYPE lg_term2 = lgammafn(term2), di_term2 = digamma(term2), tri_term2 = trigamma(term2);
            TYPE k1 = x->_a[i] + a1;
            TYPE lg_k1 = lgammafn(k1), di_k1 = x->_m1 * digamma(k1), tri_k1 = x->_m1 * x->_m1 * digamma(k1);
            TYPE k2 = x->_ta[i] - x->_a[i] + a2;
            TYPE lg_k2 = lgammafn(k2), di_k2 = m2 * digamma(k2), tri_k2 = m2 * m2 * trigamma(k2);
            f -= (lg_a0 - lg_term2 + lg_k1 - lg_a1 + lg_k2 - lg_a2);
            *g -= (di_a0 - di_term2 + di_k1 - di_a1 + di_k2 - di_a2);
            *h -= (tri_a0 - tri_term2 + tri_k1 - tri_a1 + tri_k2 - tri_a2);
        }
        return f;
    } else {
        TYPE a1 = s * x->_m1, a2 = s * m2;
        TYPE f = -(lgammafn(s) - lgammafn(a1) - lgammafn(a2)) * x->_size;
        forint(i, 0, x->_size) {
            f -= (-lgammafn(x->_ta[i] + s) + lgammafn(x->_a[i] + a1) + lgammafn(x->_ta[i] - x->_a[i] + a2));
        }
        return f;
    }
}

TYPE fval_s_inv(TYPE* g, TYPE* h, bb_t* x, TYPE s_inv) {
    TYPE m2 = 1.0 - x->_m1, s = 1.0 / s_inv;
    if (g) {
        TYPE f = 0.0; *g = 0.0; *h = 0.0;
        TYPE lg_a0 = lgammafn(s), di_a0 = digamma(s), tri_a0 = trigamma(s);
        TYPE a1 = s * x->_m1, a2 = s * m2;
        TYPE lg_a1 = lgammafn(a1), di_a1 = x->_m1 * digamma(a1), tri_a1 = x->_m1 * x->_m1 * trigamma(a1);
        TYPE lg_a2 = lgammafn(a2), di_a2 = m2 * digamma(a2), tri_a2 = m2 * m2 * trigamma(a2);
        forint(i, 0, x->_size) {
            TYPE term2 = x->_ta[i] + s;
            TYPE di_term2 = digamma(term2), tri_term2 = trigamma(term2);
            TYPE k1 = x->_a[i] + a1;
            TYPE di_k1 = x->_m1 * digamma(k1), tri_k1 = x->_m1 * x->_m1 * digamma(k1);
            TYPE k2 = x->_ta[i] - x->_a[i] + a2;
            TYPE di_k2 = m2 * digamma(k2), tri_k2 = m2 * m2 * trigamma(k2);
            f -= (lgammafn(s) - lgammafn(term2) + lgammafn(k1) - lgammafn(a1) + lgammafn(k2) - lgammafn(a2));
            TYPE tmp_g = (di_a0 - di_term2 + di_k1 - di_a1 + di_k2 - di_a2);
            TYPE s_inv_sqr = s_inv * s_inv;
            *g += s_inv_sqr * tmp_g;
            *h -= (2.0 * s_inv_sqr * s_inv * tmp_g + s_inv_sqr * s_inv_sqr *
                (tri_a0 - trigamma(term2) + tri_k1 - tri_a1 + tri_k2 - tri_a2));
        }
        return f;
    } else {
        TYPE a1 = s * x->_m1, a2 = s * m2;
        TYPE f = -(lgammafn(s) - lgammafn(a1) - lgammafn(a2)) * x->_size;
        forint(i, 0, x->_size) {
            f -= (-lgammafn(x->_ta[i] + s) + lgammafn(x->_a[i] + a1) + lgammafn(x->_ta[i] - x->_a[i] + a2));
        }
        return f;
    }
}

TYPE fval_s_equal_inv(TYPE* dx, TYPE* dxx, bb_t* x, TYPE s_inv) {
    if (dx) {
        TYPE f = 0.0; *dx = 0.0; *dxx = 0.0; TYPE _dx, _dxx;
        forint(g, 0, x->M) {
            x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g];
            x->_m1 = x->_m1_array[g];
            f += fval_s_inv(&_dx, &_dxx, x, s_inv); *dx += _dx; *dxx += _dxx;
        }
        return f;
    } else {
        TYPE f = 0.0;
        forint(g, 0, x->M) {
            x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g];
            x->_m1 = x->_m1_array[g]; f += fval_s_inv(0, 0, x, s_inv);
        }
        return f;
    }
}

void fit_m(bb_t* x, const TYPE s) {
    forint(i, 0, 200) {
        TYPE m2 = 1.0 - x->_m1, alp = s * x->_m1, bet = s * m2;
        TYPE new_m1 = -digamma(alp) * x->_size, new_m2 = -digamma(bet) * x->_size;
        forint(j, 0, x->_size) {
            new_m1 += digamma(alp + x->_a[j]);
            new_m2 += digamma(bet + x->_ta[j] - x->_a[j]);
        }
        new_m1 *= alp; new_m2 *= bet;
        TYPE tmp = new_m1 + new_m2; new_m1 /= tmp;
        if (fabs(x->_m1 - new_m1) < 1e-8) {
            x->_m1 = new_m1;
            if (x->_m1 < SMALL) x->_m1 = SMALL;
            if ((x->_m1 + SMALL) > 1.0) x->_m1 = 1.0 - SMALL;
            break;
        }
        x->_m1 = new_m1;
        if (x->_m1 < SMALL) x->_m1 = SMALL;
        if ((x->_m1 + SMALL) > 1.0) x->_m1 = 1.0 - SMALL;
    }
}

TYPE bbmle(bb_t* x, int g, TYPE alp0, TYPE bet0, TYPE* alp, TYPE* bet) {
    if (g >= 0) { x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g]; }
    else { x->_size = x->N; x->_a = x->a; x->_ta = x->ta; }
    TYPE s_inv = 1.0 / (alp0 + bet0);
    x->_m1 = alp0 * s_inv;
    forint(i, 0, 5000) {
        TYPE old_m1 = x->_m1;
        fit_m(x, 1.0 / s_inv);
        TYPE old_s_inv = s_inv, old_f = fval_s_inv(0, 0, x, old_s_inv);
        s_inv = minimize1d(fval_s_inv, x, old_s_inv, SMIN, SMAX);
        TYPE new_f = fval_s_inv(0, 0, x, s_inv);
        if (fabs(s_inv - old_s_inv) < 1e-12 && fabs(x->_m1 - old_m1) < 1e-12 && fabs(old_f - new_f) < 1e-12) break;
    }
    *alp = x->_m1 / s_inv; *bet = (1.0 - x->_m1) / s_inv;
    return -fval_s_inv(0, 0, x, s_inv);
}

TYPE bbmle_equal(bb_t* x, TYPE alp, TYPE bet) {
    TYPE alp0, bet0, s;
    if (alp < 0.0) {
        s = 0.0;
        forint(g, 0, x->M) {
            x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g];
            bb_simple_estimate_tm(x, &alp0, &bet0);
            TYPE tmp_s = alp0 + bet0; x->_m1_array[g] = alp0 / tmp_s; s += tmp_s;
        }
        s /= (TYPE)x->M;
    } else {
        s = alp + bet;
        forint(g, 0, x->M) x->_m1_array[g] = alp / s;
    }
    TYPE s_inv = 1.0 / s;
    forint(i, 0, 5000) {
        TYPE max_diff = 0.0;
        forint(g, 0, x->M) {
            x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g];
            x->_m1 = x->_m1_array[g]; TYPE old_m1 = x->_m1;
            fit_m(x, 1.0 / s_inv);
            x->_m1_array[g] = x->_m1;
            if (max_diff < fabs(x->_m1 - old_m1)) max_diff = fabs(x->_m1 - old_m1);
        }
        TYPE old_s_inv = s_inv;
        s_inv = minimize1d(fval_s_equal_inv, x, old_s_inv, SMIN, SMAX);
        if (fabs(s_inv - old_s_inv) < 1e-12 && max_diff < 1e-12) break;
    }
    return -fval_s_equal_inv(0, 0, x, s_inv);
}

void do_bb_test(bb_t* x) {
    TYPE alp, bet;
    x->_size = x->N; x->_a = x->a; x->_ta = x->ta;
    bb_simple_estimate_tm(x, &alp, &bet);
    TYPE alp0 = alp, bet0 = bet;
    TYPE f0 = bbmle(x, -1, alp0, bet0, &alp, &bet);
    alp0 = alp; bet0 = bet;
    TYPE f;
    if (x->theta_equal > 0) {
        TYPE f_init1 = bbmle_equal(x, alp0, bet0);
        int comp = x->_m1_array[0] > x->_m1_array[1];
        TYPE f_init2 = bbmle_equal(x, -1.0, -1.0);
        if (f_init1 > f_init2) { f = f_init1; x->comp = comp; }
        else { f = f_init2; x->comp = x->_m1_array[0] > x->_m1_array[1]; }
    } else {
        f = 0.0; TYPE v0 = 0, v1 = 0;
        forint(g, 0, x->M) {
            TYPE tmp = bbmle(x, g, alp0, bet0, &alp, &bet); TYPE vv = x->_m1;
            x->_size = x->g_size[g]; x->_a = x->a + x->g_ind[g]; x->_ta = x->ta + x->g_ind[g];
            bb_simple_estimate_tm(x, &alp, &bet);
            TYPE alp_tmp, bet_tmp;
            TYPE tmp2 = bbmle(x, g, alp, bet, &alp_tmp, &bet_tmp);
            if (tmp > tmp2) { f += tmp; } else { f += tmp2; vv = x->_m1; }
            if (g == 0) v0 = vv;
            if (g == 1) v1 = vv;
        }
        x->comp = v0 > v1;
    }
    x->f0 = f0; x->f = f;
}

void bbCores(int* n_procs) {
    *n_procs = 1;
#ifdef _OPENMP
    *n_procs = omp_get_num_procs();
#endif
}

void bb(int* lK, double* a, double* ta, int* lM, int* g_size, int* g_ind,
        double* mem, int* no_threads, double* pval) {
    int verbose = *no_threads;
    if (*no_threads < 0) *no_threads = -*no_threads;
    int theta_equal = ((TYPE)mem[0] > 0 ? 1 : 0);
    TYPE tail = (TYPE)(mem[1]);
    int N = 0;
    forint(i, 0, *lM) N += g_size[i];
    int block = 2 * N + (*lM);
    TYPE* works = (TYPE*)mem;
    int thres_display = 0;
#ifdef _OPENMP
    omp_set_dynamic(0); omp_set_num_threads(*no_threads);
    if (verbose > 0) printf("Using %d thread(s), rows=%d, groups=%d, samples=%d\n", *no_threads, *lK, *lM, N);
#else
    if (verbose > 0) printf("Single thread, rows=%d, groups=%d, samples=%d\n", *lK, *lM, N);
#endif
    int k = 0, stop_sig = 0;
    #pragma omp parallel for schedule(dynamic)
    forint(i, 0, *lK) {
        if (stop_sig) continue;
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        bb_t x;
        x.N = N; x.theta_equal = theta_equal; x.M = *lM; x.g_size = g_size; x.g_ind = g_ind;
        int mem_start = thread_id * block;
        x.a = works + mem_start; x.ta = works + mem_start + N; x._m1_array = works + mem_start + 2 * N;
        int ind = i * N;
        forint(j, 0, N) { x.a[j] = (TYPE)a[ind]; x.ta[j] = (TYPE)ta[ind]; ind++; }
        do_bb_test(&x);
        double g = 2.0 * ((double)x.f - (double)x.f0);
        if (tail > 0.5) {
            pval[i] = x.comp ? pnorm(-sqrt(g), 0, 1, 0, 0) : pnorm(sqrt(g), 0, 1, 0, 0);
        } else if (tail < -0.5) {
            pval[i] = x.comp ? pnorm(-sqrt(g), 0, 1, 1, 0) : pnorm(sqrt(g), 0, 1, 1, 0);
        } else {
            pval[i] = pchisq(g, (double)*lM - 1.0, 0, 0);
        }
        #pragma omp atomic
        k++;
        if (thread_id == 0 && verbose > 0 && k > thres_display) {
            printf("%d%%\n", k * 100 / (*lK)); fflush(stdout);
            thres_display = k + (*lK) / 20;
        }
    }
    if (verbose > 0) printf("Done.\n");
}
