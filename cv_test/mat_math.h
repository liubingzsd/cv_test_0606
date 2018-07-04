#ifndef _MAT_MATH_H_
#define _MAT_MATH_H_
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "mat.h"

void mat_printf(mat_s *mat);
mat_s *mat_mul(mat_s *a,mat_s *b);
mat_s *mat_sub(mat_s *a, mat_s *b);
mat_s *mat_add(mat_s *a, mat_s *b);
void mat_sub_r(mat_s *a, mat_s *b, mat_s *c);
#endif

