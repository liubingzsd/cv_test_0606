#ifndef _IMAGE_FILTER_H_
#define _IMAGE_FILTER_H_
#include "mat.h"
#include "kernel.h"
enum boundary_type {
	NONE,
	ZERO,
	REPLICATE
};

enum smooth_type {
	GAUSSIAN
};

void image_filter_1d_x(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type);
void image_filter_1d_y(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type);
void image_filter_1d_x_none(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type);
void image_filter_1d_y_none(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type);



#endif

