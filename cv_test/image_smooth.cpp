#include "mat.h"
#include "common.h"
#include "kernel.h"
#include "image_filter.h"
#include "image_smooth.h"

void image_smooth_gaussian(mat_s *src, mat_s *dst, float sig1, float sig2)
{
	kernel_s ker;
	mat_s* tmp;
	tmp = mat_create(src->rows,src->cols,MAT_F32);
	int dim = (int)MAX(3.0f, 2.0f * 3.5f * sig1 + 1.0f);
	if (dim % 2 == 0)
	{
		dim++;
	}
	write_mat_data_to_file_f32((char *)"src.dat", src);
	kernel_create_1d(&ker,dim,KERNEL_F32);
	gauss_kernel_1d(&ker,sig1);
	tmp = mat_clone_f32(src);
	image_filter_1d_x(src, tmp, &ker, NONE);
	write_mat_data_to_file_f32((char *)"tmp.dat", src);
	dst = mat_clone_f32(tmp);
	image_filter_1d_y(tmp, dst, &ker, NONE);
	write_mat_data_to_file_f32((char *)"dst.dat", dst);
	mat_free(tmp);
	kernel_free(&ker);
}

mat_s* image_smooth_gaussian_rt(mat_s *src,float sig1, float sig2)
{
	kernel_s ker;
	mat_s* tmp,*dst;
	tmp = mat_create(src->rows, src->cols, MAT_F32);
	int dim = (int)MAX(3.0f, 2.0f * 3.5f * sig1 + 1.0f);
	if (dim % 2 == 0)
	{
		dim++;
	}
	kernel_create_1d(&ker, dim, KERNEL_F32);
	gauss_kernel_1d(&ker, sig1);
	tmp = mat_clone_f32(src);
	image_filter_1d_x(src, tmp, &ker, NONE);
	dst = mat_clone_f32(tmp);
	image_filter_1d_y(tmp, dst, &ker, NONE);
	mat_free(tmp);
	kernel_free(&ker);
	return dst;
}

void test_image_smooth()
{

}