#include "mat_math.h"

mat_s *mat_mul(mat_s *a,mat_s *b)
{
	mat_s *c = NULL;
	int i, j, l, u;
	if (a->cols != b->rows)
	{
		return c;
	}
	else
	{
		c = mat_create(a->rows, b->cols, a->type);
		float *a_buf = (float *)a->buf;
		float *b_buf = (float *)b->buf;
		float *c_buf = (float *)c->buf;
		for (i = 0; i < a->rows; i++)
		{
			for (j = 0; j < b->cols; j++)
			{
				u = i * b->rows + j;
				c_buf[u] = 0.0f;
				for (l = 0; l < a->cols; l++)
				{
					c_buf[u] += a_buf[i * a->cols + l] * b_buf[l * b->cols + j];
				}
			}
		}
	}
	return c;
}

mat_s *mat_sub(mat_s *a,mat_s *b)
{
	mat_s *c = NULL;
	int i, j, idx, row_idx;
	c = mat_create(a->rows, a->cols, MAT_F32);
	float *a_buf = (float *)a->buf;
	float *b_buf = (float *)b->buf;
	float *c_buf = (float *)c->buf;
	for (i = 0; i < a->rows; i++)
	{
		row_idx = i * a->cols;
		for (j = 0; j < a->cols; j++)
		{
			idx = row_idx + j;
			c_buf[idx] = a_buf[idx] - b_buf[idx];
		}
	}
	return c;
}
void mat_sub_r(mat_s *a, mat_s *b,mat_s *c)
{
	int i, j, idx, row_idx;
	float *a_buf = (float *)a->buf;
	float *b_buf = (float *)b->buf;
	float *c_buf = (float *)c->buf;
	for (i = 0; i < a->rows; i++)
	{
		row_idx = i * a->cols;
		for (j = 0; j < a->cols; j++)
		{
			idx = row_idx + j;
			c_buf[idx] = a_buf[idx] - b_buf[idx];
		}
	}
}

mat_s *mat_add(mat_s *a, mat_s *b)
{
	mat_s *c = NULL;
	int i, j, idx, row_idx;
	c = mat_create(a->rows, a->cols, MAT_F32);
	float *a_buf = (float *)a->buf;
	float *b_buf = (float *)b->buf;
	float *c_buf = (float *)c->buf;
	for (i = 0; i < a->rows; i++)
	{
		row_idx = i * a->cols;
		for (j = 0; j < a->cols; j++)
		{
			idx = row_idx + j;
			c_buf[idx] = a_buf[idx] + b_buf[idx];
		}
	}
	return c;
}
