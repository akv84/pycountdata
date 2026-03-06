/*
 * Inverted beta-binomial test - standalone C for Python ctypes.
 * Original: Thang V. Pham, t.pham@amsterdamumc.nl
 */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "rmath_standalone.h"

#define T_LN2 0.693147180559945309417232121458
#define ZERO_REG (TYPE)0.5
#define MAX_INNER 2000
#undef SMALL
#define SMALL 1e-12
#undef ALPHA
#define ALPHA  0.1
#undef BETA
#define BETA 0.7
#define SMAX (TYPE)1e6
#define SMIN (TYPE)1e-6
#define forint(i, a, b) for (int i=(int)(a); i<(int)(b); i++)
#define TYPE double

#include "q15.h"

typedef struct _data_t {
    TYPE* a; TYPE* b; TYPE* ta; TYPE* tb;
    int N;
    const TYPE* qz; const TYPE* qlw; int Z;
    TYPE* ll; TYPE* log_one_plus_z; TYPE* log_one_minus_z; TYPE* works;
    TYPE f0, f, alp, bet;
} data_t;

TYPE fval2(TYPE* ga, TYPE* gb, TYPE* haa, TYPE* hab, TYPE* hbb, data_t* x, TYPE alp, TYPE bet) {
    TYPE* y = x->works, *w = x->works + x->Z, *da_y = x->works + 2*x->Z;
    TYPE* db_y = x->works + 3*x->Z, *tmp = x->works + 4*x->Z;
    TYPE a1 = alp - 1.0, b1 = bet - 1.0;
    TYPE c = -(alp + bet - 1.0) * T_LN2 + lgammafn(alp + bet) - lgammafn(alp) - lgammafn(bet);
    forint(i, 0, x->Z) tmp[i] = x->qlw[i] + a1 * x->log_one_plus_z[i] + b1 * x->log_one_minus_z[i] + c;
    TYPE digamma_ab, digamma_a, digamma_b, trigamma_ab, trigamma_a, trigamma_b;
    TYPE f = 0.0;
    if (ga) {
        *ga = *gb = *haa = *hab = *hbb = 0.0;
        digamma_ab = digamma(alp + bet); digamma_a = digamma(alp); digamma_b = digamma(bet);
        trigamma_ab = trigamma(alp + bet); trigamma_a = trigamma(alp); trigamma_b = trigamma(bet);
    }
    TYPE* ptr = x->ll;
    forint(n, 0, x->N) {
        TYPE max_y = tmp[0] + ptr[0];
        forint(i, 0, x->Z) { y[i] = tmp[i] + *ptr++; if (max_y < y[i]) max_y = y[i]; }
        TYPE se = 0.0;
        forint(i, 0, x->Z) se += exp(y[i] - max_y);
        TYPE fx = log(se) + max_y;
        f -= fx;
        if (ga) {
            TYPE da_f = 0.0, db_f = 0.0;
            forint(i, 0, x->Z) {
                w[i] = exp(y[i] - fx);
                da_y[i] = x->log_one_plus_z[i] - T_LN2 + digamma_ab - digamma_a;
                db_y[i] = x->log_one_minus_z[i] - T_LN2 + digamma_ab - digamma_b;
                da_f += w[i] * da_y[i]; db_f += w[i] * db_y[i];
            }
            TYPE daa_f = trigamma_ab - trigamma_a, dab_f = trigamma_ab, dbb_f = trigamma_ab - trigamma_b;
            forint(i, 0, x->Z) {
                TYPE da_w = w[i] * (da_y[i] - da_f), db_w = w[i] * (db_y[i] - db_f);
                daa_f += da_w * da_y[i]; dab_f += da_w * db_y[i]; dbb_f += db_w * db_y[i];
            }
            *ga -= da_f; *gb -= db_f; *haa -= daa_f; *hab -= dab_f; *hbb -= dbb_f;
        }
    }
    return f;
}

TYPE ibb_fval(TYPE* g, TYPE* h, data_t* x, TYPE alp, TYPE bet, int derivative) {
    TYPE* y = x->works, *w = x->works + x->Z, *d_y = x->works + 2*x->Z, *tmp = x->works + 3*x->Z;
    TYPE a1 = alp - 1.0, b1 = bet - 1.0;
    TYPE c = -(alp + bet - 1.0) * T_LN2 + lgammafn(alp + bet) - lgammafn(alp) - lgammafn(bet);
    forint(i, 0, x->Z) tmp[i] = x->qlw[i] + a1 * x->log_one_plus_z[i] + b1 * x->log_one_minus_z[i] + c;
    TYPE digamma_ab, digamma_a, digamma_b, trigamma_ab, trigamma_a, trigamma_b;
    TYPE f = 0.0;
    if (g) {
        *g = *h = 0.0;
        digamma_ab = digamma(alp + bet); digamma_a = digamma(alp); digamma_b = digamma(bet);
        trigamma_ab = trigamma(alp + bet); trigamma_a = trigamma(alp); trigamma_b = trigamma(bet);
    }
    TYPE* ptr = x->ll;
    forint(n, 0, x->N) {
        TYPE max_y = tmp[0] + ptr[0];
        forint(i, 0, x->Z) { y[i] = tmp[i] + *ptr++; if (max_y < y[i]) max_y = y[i]; }
        TYPE se = 0.0;
        forint(i, 0, x->Z) se += exp(y[i] - max_y);
        TYPE fx = log(se) + max_y;
        f -= fx;
        if (g) {
            TYPE df = 0.0;
            if (derivative == 0) {
                for (int i = 0; i < x->Z; i++) {
                    w[i] = exp(y[i] - fx);
                    d_y[i] = x->log_one_plus_z[i] - T_LN2 + digamma_ab - digamma_a;
                    df += w[i] * d_y[i];
                }
                TYPE ddf = trigamma_ab - trigamma_a;
                for (int i = 0; i < x->Z; i++) { ddf += w[i] * (d_y[i] - df) * d_y[i]; }
                *g -= df; *h -= ddf;
            } else if (derivative == 1) {
                for (int i = 0; i < x->Z; i++) {
                    w[i] = exp(y[i] - fx);
                    d_y[i] = x->log_one_minus_z[i] - T_LN2 + digamma_ab - digamma_b;
                    df += w[i] * d_y[i];
                }
                TYPE ddf = trigamma_ab - trigamma_b;
                for (int i = 0; i < x->Z; i++) { ddf += w[i] * (d_y[i] - df) * d_y[i]; }
                *g -= df; *h -= ddf;
            } else if (derivative == 2) {
                for (int i = 0; i < x->Z; i++) {
                    w[i] = exp(y[i] - fx);
                    d_y[i] = x->log_one_plus_z[i] + x->log_one_minus_z[i] - 2.0 * T_LN2
                             + 2.0 * digamma_ab - digamma_a - digamma_b;
                    df += w[i] * d_y[i];
                }
                TYPE ddf = 4 * trigamma_ab - trigamma_a - trigamma_b;
                for (int i = 0; i < x->Z; i++) { ddf += w[i] * (d_y[i] - df) * d_y[i]; }
                *g -= df; *h -= ddf;
            }
        }
    }
    return f;
}

TYPE fval_ab(TYPE* g, TYPE* h, data_t* x, TYPE p) {
    if (g) return ibb_fval(g, h, x, p, p, 2);
    else return ibb_fval(g, h, x, p, p, -1);
}

TYPE ibb_find_eta(TYPE* v, TYPE g, TYPE h, TYPE p, TYPE lb, TYPE ub) {
    if (fabs(h) < SMALL) h = (h > 0 ? SMALL : -SMALL);
    *v = -g / h;
    if (g * (*v) > 0) *v = -(*v);
    TYPE eta = 1.0;
    if (p + *v >= ub) eta = (ub - SMALL - p) / *v;
    if (p + *v <= lb) eta = (lb + SMALL - p) / *v;
    if (eta < 0.0) eta = 0.0;
    return eta;
}

void update_ab(TYPE* new_f, data_t* x, TYPE* a, TYPE* b, TYPE f, TYPE fcond,
               TYPE ga, TYPE gb, TYPE eta, TYPE va, TYPE vb,
               TYPE ua, TYPE la, TYPE ub, TYPE lb, TYPE ab_bound) {
    TYPE lambda = ga * va + gb * vb;
    *new_f = fval2(0, 0, 0, 0, 0, x, *a + eta * va, *b + eta * vb) * fcond;
    while ((*new_f - f) > ALPHA * lambda * eta) {
        eta *= BETA;
        *new_f = fval2(0, 0, 0, 0, 0, x, *a + eta * va, *b + eta * vb) * fcond;
        if (eta < 1e-15) { eta = 0.0; *new_f = f; break; }
    }
    *a += eta * va; *b += eta * vb;
}

void nr2b_projection(TYPE(*fval2_fn)(TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, data_t*, TYPE, TYPE),
    data_t* x, TYPE* a, TYPE la, TYPE ua, TYPE* b, TYPE lb, TYPE ub,
    TYPE ab_bound, TYPE fcond) {
    TYPE old_a = *a, old_b = *b;
    int cc = 0, slide_ok = 1, newton_ok = 1;
    forint(iter, 0, MAX_INNER) {
        TYPE ga, gb, haa, hab, hbb;
        TYPE f = fval2_fn(&ga, &gb, &haa, &hab, &hbb, x, *a, *b);
        f *= fcond; ga *= fcond; gb *= fcond; haa *= fcond; hab *= fcond; hbb *= fcond;
        if (fabs(ga) < 1e-20 && fabs(gb) < 1e-20) break;
        TYPE va = 0.0, vb = 0.0, new_f = f;
        TYPE det = haa * hbb - hab * hab, eta = 1.0;
        TYPE da = gb - ga, db = ga - gb;
        TYPE d_ddf_d = da * (haa * da + hab * db) + db * (hab * da + hbb * db);
        TYPE df_d = ga * da + gb * db;
        if ((*a + *b - ab_bound) < 1e-9 && d_ddf_d > 1e-8 && slide_ok) {
            TYPE delta = -df_d / d_ddf_d; if (delta < 0) delta = -delta;
            va = da * delta; vb = db * delta;
            if (*a + eta * va >= ua) eta = (ua - SMALL - *a) / va;
            if (*a + eta * va <= la) eta = (la + SMALL - *a) / va;
            if (*b + eta * vb >= ub) eta = (ub - SMALL - *b) / vb;
            if (*b + eta * vb <= lb) eta = (lb + SMALL - *b) / vb;
            update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);
            if (fabs(new_f - f) < 1e-10 && fabs(*a - old_a) < 1e-10 && fabs(*b - old_b) < 1e-10) slide_ok = 0;
            else { slide_ok = 1; newton_ok = 1; }
            old_a = *a; old_b = *b; continue;
        }
        if (haa > 1e-30 && det > 1e-30 && newton_ok) {
            va = -(hbb * ga - hab * gb) / det; vb = -(-hab * ga + haa * gb) / det;
            if (*a + eta * va >= ua) eta = (ua - SMALL - *a) / va;
            if (*a + eta * va <= la) eta = (la + SMALL - *a) / va;
            if (*b + eta * vb >= ub) eta = (ub - SMALL - *b) / vb;
            if (*b + eta * vb <= lb) eta = (lb + SMALL - *b) / vb;
            if (*a + eta * va + *b + eta * vb <= ab_bound) eta = (ab_bound + SMALL - *a - *b) / (va + vb);
            update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);
            if (fabs(new_f - f) < 1e-10 && fabs(*a - old_a) < 1e-10 && fabs(*b - old_b) < 1e-10) newton_ok = 0;
            else { slide_ok = 1; newton_ok = 1; }
            old_a = *a; old_b = *b; continue;
        }
        TYPE eta_a = ibb_find_eta(&va, ga, haa, *a, la > (ab_bound - *b) ? la : (ab_bound - *b), ua);
        TYPE eta_b = ibb_find_eta(&vb, gb, hbb, *b, lb > (ab_bound - *a) ? lb : (ab_bound - *a), ub);
        if (fabs(eta_a * va) > fabs(eta_b * vb)) { eta = eta_a; vb = 0.0; }
        else { eta = eta_b; va = 0.0; }
        update_ab(&new_f, x, a, b, f, fcond, ga, gb, eta, va, vb, ua, la, ub, lb, ab_bound);
        if (fabs(new_f - f) < 1e-10 && fabs(*a - old_a) < 1e-10 && fabs(*b - old_b) < 1e-10) cc++;
        else { cc = 0; slide_ok = 1; newton_ok = 1; }
        if (cc > 1) break;
        old_a = *a; old_b = *b;
    }
}

TYPE nr(TYPE(*fval_fn)(TYPE*, TYPE*, data_t*, TYPE), data_t* x, TYPE p0, TYPE lb, TYPE ub) {
    TYPE p = p0, old_p = p;
    int cc = 0;
    for (int i = 0; i < MAX_INNER; i++) {
        TYPE f, g, h;
        f = fval_fn(&g, &h, x, p);
        TYPE v, eta = ibb_find_eta(&v, g, h, p, lb, ub);
        TYPE lambda = g * v, fnew = fval_fn(0, 0, x, p + eta * v);
        while ((fnew - f) > ALPHA * lambda * eta) {
            eta *= BETA; fnew = fval_fn(0, 0, x, p + eta * v);
            if (eta < 1e-15) { eta = 0.0; fnew = f; break; }
        }
        p += eta * v;
        if (fabs(fnew - f) < 1e-10 && fabs(p - old_p) < 1e-10) cc++; else cc = 0;
        if (cc > 1) break;
        old_p = p;
    }
    return p;
}

void do_ibb_test(data_t* x) {
    TYPE ab_bound = x->works[0];
    TYPE lower_bound = 1.0, upper_bound = 1e5;
    TYPE tmp = (lower_bound > (ab_bound / 2.0) ? lower_bound : (ab_bound / 2.0));
    TYPE ab0 = tmp + 1.0; ab0 = nr(fval_ab, x, ab0, tmp, upper_bound);
    TYPE fab0 = ibb_fval(0, 0, x, ab0, ab0, -1);
    TYPE ab1 = upper_bound - 1.0; ab1 = nr(fval_ab, x, ab1, tmp, upper_bound);
    TYPE fab1 = ibb_fval(0, 0, x, ab1, ab1, -1);
    TYPE ab2 = (upper_bound + tmp) / 2.0; ab2 = nr(fval_ab, x, ab2, tmp, upper_bound);
    TYPE fab2 = ibb_fval(0, 0, x, ab2, ab2, -1);
    TYPE f0 = fab0, ab = ab0;
    if (f0 > fab1) { f0 = fab1; ab = ab1; }
    if (f0 > fab2) { f0 = fab2; ab = ab2; }
    TYPE alp = ab, bet = ab;
    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, 1.0);
    TYPE f = fval2(0, 0, 0, 0, 0, x, alp, bet);
    TYPE fmin = f, amin = alp, bmin = bet;
    alp = tmp + 1.0; bet = upper_bound - 1.0;
    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, 1.0);
    f = fval2(0, 0, 0, 0, 0, x, alp, bet);
    if (f < fmin) { fmin = f; amin = alp; bmin = bet; }
    alp = upper_bound - 1.0; bet = tmp + 1.0;
    nr2b_projection(fval2, x, &alp, lower_bound, upper_bound, &bet, lower_bound, upper_bound, ab_bound, 1.0);
    f = fval2(0, 0, 0, 0, 0, x, alp, bet);
    if (f < fmin) { fmin = f; amin = alp; bmin = bet; }
    x->f0 = -f0; x->f = -fmin; x->alp = amin; x->bet = bmin;
}

void ibb(int* lK, double* aa, double* bb_arr, double* taa, double* tbb,
         int* lN, double* mem, int* no_threads, double* pval, double* fc) {
    int verbose = *no_threads;
    if (*no_threads < 0) *no_threads = -*no_threads;
    TYPE lower_bound = (TYPE)(mem[0]), tail = (TYPE)(mem[1]);
    int init_block = 3 * _Z;
    int block = 4 * (*lN) + (*lN) * _Z + 5 * _Z;
    TYPE* works = (TYPE*)mem;
    TYPE* log_one_minus_z = works, *log_one_plus_z = works + _Z, *phi = works + 2 * _Z;
    int thres_display = 0;
    forint(i, 0, _Z) {
        log_one_minus_z[i] = log1p(-_qz[i]);
        log_one_plus_z[i] = log1p(_qz[i]);
        phi[i] = (1.0 + _qz[i]) / (1.0 - _qz[i]);
    }
#ifdef _OPENMP
    omp_set_dynamic(0); omp_set_num_threads(*no_threads);
    if (verbose > 0) printf("Using %d thread(s), rows=%d, pairs=%d\n", *no_threads, *lK, *lN);
#else
    if (verbose > 0) printf("Single thread, rows=%d, pairs=%d\n", *lK, *lN);
#endif
    int k = 0, stop_sig = 0;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < *lK; i++) {
        if (stop_sig) continue;
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        data_t x;
        x.N = *lN; x.qz = _qz; x.qlw = _qw; x.Z = _Z;
        x.log_one_minus_z = log_one_minus_z; x.log_one_plus_z = log_one_plus_z;
        int mem_start = init_block + thread_id * block;
        x.a = works + mem_start; x.b = works + mem_start + (*lN);
        x.ta = works + mem_start + 2 * (*lN); x.tb = works + mem_start + 3 * (*lN);
        x.ll = works + mem_start + 4 * (*lN);
        x.works = works + mem_start + 4 * (*lN) + (*lN) * _Z;
        int ind = i * (*lN);
        TYPE* ptr = x.ll;
        forint(j, 0, (*lN)) {
            x.a[j] = (TYPE)aa[ind]; x.b[j] = (TYPE)bb_arr[ind];
            x.ta[j] = (TYPE)taa[ind]; x.tb[j] = (TYPE)tbb[ind];
            if ((x.a[j] + x.b[j]) > 0) {
                TYPE factor = x.ta[j] / x.tb[j];
                forint(kk, 0, _Z) { TYPE p = phi[kk] / (phi[kk] + factor); *ptr++ = x.a[j] * log1p(-p) + x.b[j] * log(p); }
            } else {
                forint(kk, 0, _Z) { TYPE p = phi[kk] / (phi[kk] + 1.0); *ptr++ = ZERO_REG * (log1p(-p) + log(p)); }
            }
            ind++;
        }
        x.works[0] = lower_bound;
        do_ibb_test(&x);
        fc[i] = (double)(x.alp / x.bet);
        double g = 2.0 * ((double)x.f - (double)x.f0);
        if (tail > 0.5) {
            pval[i] = x.alp > x.bet ? pnorm(sqrt(g), 0, 1, 0, 0) : pnorm(-sqrt(g), 0, 1, 0, 0);
        } else if (tail < -0.5) {
            pval[i] = x.alp > x.bet ? pnorm(sqrt(g), 0, 1, 1, 0) : pnorm(-sqrt(g), 0, 1, 1, 0);
        } else {
            pval[i] = pchisq(g, 1.0, 0, 0);
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
