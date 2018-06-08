#include "mat.h"
#include "common.h"
#include "kernel.h"
#include "mat_math.h"

#define GAUSS_KERN  3.5


void mat_normalize(mat_s *mat)
{
	float *mat_buf = (float *)mat->buf;
	uint8_t i, j;
	float sum = 0.0f;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			sum += mat_buf[i * mat->cols + j];
		}
	}
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			mat_buf[i * mat->cols + j] /= sum;
		}
	}

}

mat_s *gaussian_kernel_2d(float sigma)
{
	int dim = (int)max(3.0f, 2.0f * GAUSS_KERN * sigma + 1.0f);
	if (dim % 2 == 0)
	{
		dim++;
	}
	mat_s *mat = mat_create(dim, dim, MAT_F32);
	float *mat_buf = (float *)mat->buf;
	float s2 = sigma * sigma;
	int c = dim / 2;
	float m = (float)(1.0f / (sqrt(2.0f * PI) * sigma));
	for (int i = 0; i < (dim + 1) / 2; i++)
	{
		for (int j = 0; j < (dim + 1) / 2; j++)
		{
			float v = (float)(m * exp(-(1.0f * i * i + 1.0f * j * j) / (2.0f * s2)));
			mat_buf[(c + i)*mat->cols + c + j] = v;
			mat_buf[(c - i)*mat->cols + c + j] = v;
			mat_buf[(c + i)*mat->cols + c - j] = v;
			mat_buf[(c - i)*mat->cols + c - j] = v;
		}
	}
	mat_normalize(mat);
	return mat;
}


void kernel_create(kernel_s *kernel, uint8_t size, enum kernel_data_type type)
{
	kernel->size = size;
	kernel->scale = 1;
	if (type == KERNEL_U8 || type == KERNEL_S8)
	{
		kernel->buf_size = sizeof(uint8_t) * size * size;
	}
	else if (type == KERNEL_F32)
	{
		kernel->buf_size = sizeof(float) * size * size;
	}
	kernel->buf = malloc(kernel->buf_size);
	if (kernel->buf)
	{
		//		printf("kernel buf allocated at 0x%x \n", (unsigned int)kernel->buf);
	}
	else
	{
		printf("___func___ kernel_create \n");
		printf("kernel buf allocated failed \n");
	}
}


void kernel_create_1d(kernel_s *kernel, uint8_t size, enum kernel_data_type type)
{
	kernel->size = size;
	kernel->scale = 1;
	if (type == KERNEL_U8 || type == KERNEL_S8)
	{
		kernel->buf_size = sizeof(uint8_t) * size;
	}
	else if (type == KERNEL_F32)
	{
		kernel->buf_size = sizeof(float) * size;
	}
	kernel->buf = malloc(kernel->buf_size);
	if (kernel->buf)
	{
		//		printf("kernel buf allocated at 0x%x \n", (unsigned int)kernel->buf);
	}
	else
	{
		printf("___func___ kernel_create \n");
		printf("kernel buf allocated failed \n");
	}
}

void gauss_kernel_1d(kernel_s *kernel, float sigma)
{
	uint8_t i, x_cen, xtmp;
	float sig, sum = 0, sum1 = 0, tmp;
	float *kernel_buf = (float *)kernel->buf;
	x_cen = (kernel->size - 1) / 2;
	sig = 2 * sigma * sigma;
	for (i = 0; i < kernel->size; i++)
	{
		xtmp = (i - x_cen) * (i - x_cen);
		tmp = *(kernel_buf + i) = (float)exp(-xtmp / sig) / (float)(sigma * sqrt(2 * PI));
		sum += tmp;

	}
	/*normalized to 1*/
	for (i = 0; i < kernel->size; i++)
	{

		*(kernel_buf + i) /= sum;
	}
}

void gauss_kernel_2d(kernel_s *kernel, float sigma1, float sigma2)
{
	uint8_t x_cen, y_cen, i, j, xtmp, ytmp;
	float sig, sum = 0, sum1 = 0, tmp;
	float *kernel_buf = (float *)kernel->buf;
	x_cen = y_cen = (kernel->size - 1) / 2;
	sig = 2 * sigma1 * sigma2;
	for (i = 0; i < kernel->size; i++)
	{
		xtmp = (i - x_cen) * (i - x_cen);
		for (j = 0; j < kernel->size; j++)
		{
			ytmp = (j - y_cen) * (j - y_cen);
			tmp = *(kernel_buf + i * kernel->size + j) = (float)exp(-(xtmp + ytmp) / sig) / (sig);
			sum += tmp;
		}
	}
	/*normalized to 1*/
	for (i = 0; i < kernel->size; i++)
	{
		for (j = 0; j < kernel->size; j++)
		{
			*(kernel_buf + i * kernel->size + j) /= sum;
		}
	}
}

void kernel_free(kernel_s *kernel)
{
	if (kernel->buf != NULL)
	{
		free(kernel->buf);
	}
	kernel->buf = NULL;
}

void test_filter_kernel()
{
#define  M  7
	kernel_s kernel, kernel_1d;
	float sigma1 = 2;
	float sigma2 = 2;
	float gauss[M * M] = { 0 };
	printf("kernel_1d \n");
	kernel_create_1d(&kernel_1d, M, KERNEL_F32);
	gauss_kernel_1d(&kernel_1d, sigma1);
	mat_printf_f32((float *)kernel_1d.buf, 1, M);
	printf("kernel_2d \n");
	kernel_create(&kernel, M, KERNEL_F32);
	gauss_kernel_2d(&kernel, sigma1, sigma2);
	mat_printf_f32((float *)kernel.buf, M, M);
	kernel_free(&kernel);

}