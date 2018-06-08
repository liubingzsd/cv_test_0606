#include "mat.h"
#include "mat_math.h"
#include "file_operation.h"

void mat_printf_u8(uint8_t *mat, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%d ", *(mat + i * n + j));
		}
		printf("\n");
	}
}

void mat_printf_u16(uint16_t *mat, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%d ", *(mat + i * n + j));
		}
		printf("\n");
	}
}

void mat_printf_u32(uint32_t *mat, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%d ", *(mat + i * n + j));
		}
		printf("\n");
	}
}

void mat_printf_f32(float *mat, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%f ", *(mat + i * n + j));
		}
		printf("\n");
	}
}

void mat_printf(mat_s *mat)
{
	if (mat->type == MAT_U8 || mat->type == MAT_S8)
	{
		mat_printf_u8((uint8_t *)mat->buf, mat->rows, mat->cols);
	}
	else if (mat->type == MAT_U16 || mat->type == MAT_S16)
	{
		mat_printf_u16((uint16_t *)mat->buf, mat->rows, mat->cols);
	}
	else if (mat->type == MAT_U32 || mat->type == MAT_S32)
	{
		mat_printf_u32((uint32_t *)mat->buf, mat->rows, mat->cols);
	}
	else
	{
		mat_printf_f32((float *)mat->buf, mat->rows, mat->cols);
	}
}


mat_s *mat_create(uint16_t rows, uint16_t cols, enum mat_type type)
{
	mat_s *mat  = (mat_s *)malloc(sizeof(mat_s));
	mat->type = type;
	mat->rows = rows;
	mat->cols = cols;

	if (type == MAT_S8 || type == MAT_U8)
	{
		mat->buf_size = sizeof(uint8_t) * rows * cols;
	}
	else if (type == MAT_U32 || type == MAT_S32 || type == MAT_F32)
	{
		mat->buf_size = sizeof(float) * rows * cols;
	}
	else
	{
		mat->buf_size = sizeof(int16_t) * rows * cols;
	}

	mat->buf = malloc(mat->buf_size);
	if (mat->buf)
	{
	}
	else
	{
		printf("___func___ image_create \n");
		printf("mat buf allocated failed \n");
	}
	return mat;
}

void mat_free(mat_s *mat)
{
	if (mat->buf != NULL)
	{
		free(mat->buf);
		mat->buf = NULL;
	}
}

void mat_clone_f32(mat_s *src,mat_s *dst)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = x_idx + j;
			dst_buf[idx] = src_buf[idx];
		}
	}
}

mat_s *mat_clone_f32_test(mat_s *src)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	mat_s *dst;
	dst = mat_create(src->rows, src->cols, src->type);
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = x_idx + j;
			dst_buf[idx] = src_buf[idx];
		}
	}
	return dst;
}

mat_s *mat_clone_u8_test(mat_s *src)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	mat_s *dst;
	dst = mat_create(src->rows, src->cols, src->type);
	uint8_t *src_buf = (uint8_t *)src->buf;
	uint8_t *dst_buf = (uint8_t *)dst->buf;
	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = x_idx + j;
			dst_buf[idx] = src_buf[idx];
		}
	}
	return dst;
}

mat_s *mat_convert_f32_scale(mat_s *src, float scale)
{
	float *src_buf = (float *)src->buf;
	uint16_t i, j;
	uint32_t x_idx, idx;

	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = x_idx + j;
			src_buf[idx] = src_buf[idx] * scale;
		}
	}
	return src;
}
void mat_set_zero_f32(mat_s *mat)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	float *mat_buf = NULL;
	mat_buf = (float *)mat->buf;
	if (mat->type == MAT_F32)
	{
		for (i = 0; i < mat->rows; i++)
		{
			x_idx = i * mat->cols;
			for (j = 0; j < mat->cols; j++)
			{
				idx = x_idx + j;
				mat_buf[idx] = 0.0f;
			}
		}
	}

}


void mat_set_val_int(mat_s *mat, int val)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	if (mat->type == MAT_U8 || mat->type == MAT_S8)
	{
		uint8_t *mat_buf = (uint8_t *)mat->buf;
		for (i = 0; i < mat->rows; i++)
		{
			x_idx = i * mat->cols;
			for (j = 0; j < mat->cols; j++)
			{
				idx = x_idx + j;
				mat_buf[idx] = val;
			}
		}
	}
	else if (mat->type == MAT_U16 || mat->type == MAT_S16)
	{
		uint16_t *mat_buf = (uint16_t *)mat->buf;
		for (i = 0; i < mat->rows; i++)
		{
			x_idx = i * mat->cols;
			for (j = 0; j < mat->cols; j++)
			{
				idx = x_idx + j;
				mat_buf[idx] = val;
			}
		}
	}
	else
	{
		uint32_t *mat_buf = (uint32_t *)mat->buf;
		for (i = 0; i < mat->rows; i++)
		{
			x_idx = i * mat->cols;
			for (j = 0; j < mat->cols; j++)
			{
				idx = x_idx + j;
				mat_buf[idx] = val;
			}
		}
	}
}

void mat_set_val_f32(mat_s *mat, float val)
{
	uint16_t i, j;
	uint32_t idx, x_idx;
	float *mat_buf = (float *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		x_idx = i * mat->cols;
		for (j = 0; j < mat->cols; j++)
		{
			idx = x_idx + j;
			mat_buf[idx] = val;
		}
	}
}

void mat_resize_double_f32(mat_s *src, mat_s *dst)
{
	uint16_t i, j;
	uint32_t x_idx, idx;
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	if (dst->rows != 2 * src->rows || dst->cols != 2 * src->cols)
	{
		printf("reconfirm the size of dst_image and src_image \n");
	}
	else
	{
		for (i = 0; i < dst->rows; i++)
		{
			x_idx = i * dst->cols;
			for (j = 0; j < dst->cols; j++)
			{
				idx = x_idx + j;
				dst_buf[idx] = src_buf[(i / 2)*src->cols + (j / 2)];
			}
		}
		for (i = 0; i < dst->rows; i += 2)
		{
			x_idx = i * dst->cols;
			for (j = 1; j < dst->cols - 1; j += 2)
			{
				idx = x_idx + j;
				dst_buf[idx] = 0.5f * (src_buf[(i / 2)*src->cols + (j / 2)] + (src_buf[(i / 2)*src->cols + (j / 2) + 1]));
			}
		}

		for (i = 1; i < dst->rows - 1; i += 2)
		{
			x_idx = i * dst->cols;
			for (j = 0; j < dst->cols; j += 2)
			{
				idx = x_idx + j;
				dst_buf[idx] = 0.5f * (src_buf[(i / 2)*src->cols + (j / 2)] + (src_buf[((i / 2) + 1)*src->cols + (j / 2)]));
			}
		}

		for (i = 1; i < dst->rows - 1; i += 2)
		{
			x_idx = i * dst->cols;
			for (j = 0; j < dst->cols - 1; j += 2)
			{
				idx = x_idx + j;
				dst_buf[idx] = 0.25f * (src_buf[(i / 2)*src->cols + (j / 2)] + (src_buf[(i / 2) *src->cols + (j / 2)]));
			}
		}
	}
}

mat_s *mat_resize_half_f32(mat_s *src)
{
	uint16_t i, j;
	uint16_t row, col;
	row = src->rows / 2;
	col = src->cols / 2;
	uint32_t idx, x_idx;
	mat_s *img_new = mat_create(row, col, MAT_F32);
	float *img_new_buf = (float *)img_new->buf;
	float *src_buf = (float *)src->buf;
	for (i = 0; i < row; i++)
	{
		x_idx = i * img_new->cols;
		for (j = 0; j < col; j++)
		{
			idx = x_idx + j;
			img_new_buf[idx] = src_buf[(2 * i)*src->cols + (2 * j)];
		}
	}
	return img_new;
}

void mat_convert_u8_to_f32(mat_s *src, mat_s *dst)
{
	uint8_t *src_buf = (uint8_t *)src->buf;
	float *dst_buf = (float *)dst->buf;
	uint16_t i, j;
	uint32_t x_idx, idx;
//	write_mat_data_to_file_u8("src_u8.dat", src);

	for (i = 0; i < src->rows; i++)
	{
		x_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = x_idx + j;
			dst_buf[idx] = src_buf[idx];
		}
	}
//	write_mat_data_to_file_f32("f32.dat",dst);
}


/*
mat_s *mat_convert_u8_to_f32(mat_s *src)
{
	mat_s *dst = NULL;
	uint16_t i, j;
	uint32_t row_idx, idx;
	dst = mat_create(src->rows, src->cols, MAT_F32);
	uint8_t *src_buf = (uint8_t *)src->buf;
	float *dst_buf = (float *)dst->buf;
	for (i = 0; i < src->rows; i++)
	{
		row_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = row_idx + j;
			dst_buf[row_idx] = (float)src_buf[idx];
		}
	}
	write_mat_data_to_file_f32("test.dat",dst);
	return dst;
}
*/

void mat_convert_f32_to_u8(mat_s *src,mat_s *dst)
{
	uint8_t *dst_buf = (uint8_t *)dst->buf;
	float *src_buf = (float *)src->buf;
	uint16_t i, j;
	uint32_t row_idx, idx;
//	write_mat_data_to_file_f32("f32_in.dat", src);
	for (i = 0; i < src->rows; i++)
	{
		row_idx = i * src->cols;
		for (j = 0; j < src->cols; j++)
		{
			idx = row_idx + j;
			dst_buf[idx] = src_buf[idx];
		}
	}
//	write_mat_data_to_file_u8("u8.dat", dst);
}

/*
mat_s *mat_convert_f32_to_u8(mat_s *src)
{
	mat_s *dst = NULL;
	dst = mat_create(src->rows,src->cols,MAT_U8);
	uint8_t *dst_buf = (uint8_t *)dst->buf;
	float *src_buf = (float *)src->buf;
	uint16_t i, j;
	uint32_t row_idx,idx;
	for (i = 0; i < src->rows; i++)
	{
		row_idx = i * src->rows;
		for (j = 0; j < src->cols; j++)
		{
			idx = row_idx + j;
			dst_buf[idx] = (uint8_t)(src_buf[idx]);
		}
	}
	return dst;
}
*/
/*
void mat_convert_f32_scale(mat_s *src, mat_s *dst, float scale)
{
	float *src_buf = (float *)src->buf;
	float *dst_buf = (float *)dst->buf;
	uint16_t i, j;
	uint32_t x_idx, idx;
	if (src->rows != dst->rows || src->cols != dst->cols)
	{
		printf("src_image is not equal to dst_image \n");
	}
	else
	{
		for (i = 0; i < dst->rows; i++)
		{
			x_idx = i * dst->cols;
			for (j = 0; j < dst->cols; j++)
			{
				idx = x_idx + j;
				dst_buf[idx] = src_buf[idx] * scale;
			}
		}
	}
}
*/

void write_mat_data_to_file_u8(char *str,mat_s *mat)
{
	FILE *fp = fopen(str, "w+");
	if (fp == NULL)
	{
		printf("NO file \n");
	}
	else
	{
		write_data_to_file_u8(fp, mat->buf, mat->rows, mat->cols);
	}
	fclose(fp);
}

void write_mat_data_to_file_f32(char *str, mat_s *mat)
{
	FILE *fp = fopen(str, "w+");
	if (fp == NULL)
	{
		printf("NO file \n");
	}
	else
	{
		write_data_to_file_f32(fp, mat->buf, mat->rows, mat->cols);
	}
	fclose(fp);
}


mat_s *mat_from_data_u8(uint8_t *data,int rows,int cols)
{
	mat_s *mat = NULL;
	int i, j,idx;
	mat = mat_create(rows,cols,MAT_U8);
	uint8_t *mat_buf = (uint8_t *)mat->buf;
	for (i = 0; i < mat->rows; i++)
	{
		for (j = 0; j < mat->cols; j++)
		{
			idx = i * mat->cols + j;
			mat_buf[idx] = data[idx];
		}
	}
	return mat;
}

void test_mat()
{
	int m = 10;
	int n = 10;
	mat_s *mat_u8 = NULL, *mat_u16 = NULL, *mat_u32 = NULL, *mat_f32 = NULL;
	mat_f32 = mat_create(m, n, MAT_F32);
	mat_set_val_f32(mat_f32, 7890.5643f);
	mat_printf(mat_f32);
	mat_free(mat_f32);
	mat_u32 = mat_create(m, n, MAT_U32);
	mat_set_val_int(mat_u32, 7890);
	mat_printf(mat_u32);
	mat_free(mat_u32);
}
