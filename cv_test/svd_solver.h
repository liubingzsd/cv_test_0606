#ifndef _SVD_TEST_H_
#define _SVD_TEST_H_
#include "mat.h"
void test_svd_inv();
void test_svd_solver();
mat_s* svd_inv_f32(mat_s* mat);
mat_s* svd_inv_f64(mat_s* mat);
void mat_svd_decmp(mat_s* mat, double *mu, double* mv);
bool mat_svd_solver_f64(mat_s* a, double *b, double *x);
#endif
