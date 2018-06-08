#ifndef _KERNEL_H_
#define _KERNEL_H_
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "mat.h"
enum kernel_data_type {
	KERNEL_U8,
	KERNEL_S8,
	KERNEL_F32
};
enum kernel_dim_type {
	KERNEL_1D,
	KERNEL_2D,
};
typedef struct
{
	enum kernel_data_type type;
	uint8_t size;
	float scale;
	uint32_t buf_size;
	void *buf;
}kernel_s;

void kernel_create_1d(kernel_s *kernel, uint8_t size, enum kernel_data_type type);
void kernel_create(kernel_s *kernel, uint8_t size, enum kernel_data_type type);
void gauss_kernel_2d(kernel_s *kernel, float sigma1, float sigma2);
void gauss_kernel_1d(kernel_s *kernel, float sigma);
void kernel_free(kernel_s *kernel);
void test_filter_kernel();
mat_s *gaussian_kernel_2d(float sigma);

#endif