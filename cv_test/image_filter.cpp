#include "mat.h"
#include "mat_math.h"
#include "image_filter.h"

void image_filter_1d_x(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type)
{
	uint16_t i, j;
	uint8_t j_ker, border;
	uint32_t *img_idx = NULL, x_idx;
	float tmp;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	float *kernel_buf = (float *)kernel->buf;
	img_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	border = (kernel->size - 1) >> 1;
	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = border; j < src->cols - border; j++)
		{
			tmp = 0;
			for (j_ker = 0; j_ker < kernel->size; j_ker++)
			{
				img_idx[j_ker] = x_idx + j + j_ker - border;
			}
			for (j_ker = 0; j_ker < kernel->size; j_ker++)
			{
				tmp += kernel_buf[j_ker] * src_buf[img_idx[j_ker]];
			}
			dst_buf[i * dst->cols + j] = tmp;
//			printf("tmp = %f dst_buf = %f \n", tmp, dst_buf[i * dst->cols + j]);
		}
	}
}

void image_filter_1d_y(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type)
{
	uint16_t i, j;
	uint8_t i_ker, border;
	uint32_t *img_idx = NULL, *x_idx = NULL;
	float tmp;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	float *kernel_buf = (float *)kernel->buf;
	img_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	x_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	border = (kernel->size - 1) >> 1;
	for (i = border; i < src->rows - border; i++)
	{
		for (i_ker = 0; i_ker < kernel->size; i_ker++)
		{
			x_idx[i_ker] = (i_ker + i - border) * src->cols;
		}

		for (j = 0; j < src->cols; j++)
		{
			tmp = 0;
			for (i_ker = 0; i_ker < kernel->size; i_ker++)
			{
				img_idx[i_ker] = x_idx[i_ker] + j;
			}
			for (i_ker = 0; i_ker < kernel->size; i_ker++)
			{
				tmp += kernel_buf[i_ker] * src_buf[img_idx[i_ker]];
			}
			dst_buf[i * dst->cols + j] = tmp;
		}
	}
}

void image_filter_1d_x_none(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type)
{
	uint16_t i, j;
	uint8_t j_ker, border;
	uint32_t *img_idx = NULL, x_idx;
	float tmp;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	float *kernel_buf = (float *)kernel->buf;
	dst = mat_clone_f32(src);
	img_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	border = (kernel->size - 1) >> 1;
	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = border; j < src->cols - border; j++)
		{
			tmp = 0;
			for (j_ker = 0; j_ker < kernel->size; j_ker++)
			{
				img_idx[j_ker] = x_idx + j + j_ker - border;
			}
			for (j_ker = 0; j_ker < kernel->size; j_ker++)
			{
				tmp += kernel_buf[j_ker] * src_buf[img_idx[j_ker]];
			}
			dst_buf[i * dst->cols + j] = tmp;
		}
	}
}

void image_filter_1d_y_none(mat_s *src, mat_s *dst, kernel_s *kernel, enum boundary_type type)
{
	uint16_t i, j;
	uint8_t i_ker, border;
	uint32_t *img_idx = NULL, *x_idx = NULL;
	float tmp;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	float *kernel_buf = (float *)kernel->buf;
	dst = mat_clone_f32(src);
	img_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	x_idx = (uint32_t *)malloc(sizeof(uint32_t) * kernel->size);
	border = (kernel->size - 1) >> 1;
	for (i = border; i < src->rows - border; i++)
	{
		for (i_ker = 0; i_ker < kernel->size; i_ker++)
		{
			x_idx[i_ker] = (i_ker + i - border) * src->cols;
		}

		for (j = 0; j < src->cols; j++)
		{
			tmp = 0;
			for (i_ker = 0; i_ker < kernel->size; i_ker++)
			{
				img_idx[i_ker] = x_idx[i_ker] + j;
			}
			for (i_ker = 0; i_ker < kernel->size; i_ker++)
			{
				tmp += kernel_buf[i_ker] * src_buf[img_idx[i_ker]];
			}
			dst_buf[i * dst->cols + j] = tmp;
		}
	}
}


