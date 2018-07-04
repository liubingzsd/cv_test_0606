#ifndef _MAT_H_
#define _MAT_H_

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>

enum mat_type 
{
	MAT_F64,
	MAT_F32,
	MAT_U8,
	MAT_S8,
	MAT_U16,
	MAT_S16,
	MAT_U32,
	MAT_S32
};

typedef struct
{
	enum mat_type type;
	uint16_t rows;
	uint16_t cols;
	uint32_t buf_size;
	void *buf;
}mat_s;


typedef struct
{
	enum mat_type type;
	uint16_t rows;
	uint16_t cols;
	uint32_t buf_size;
	uint8_t *buf;
}mat_u8_s;

typedef struct
{
	enum mat_type type;
	uint16_t rows;
	uint16_t cols;
	uint32_t buf_size;
	uint16_t *buf;
}mat_u16_s;

typedef struct
{
	enum mat_type type;
	uint16_t rows;
	uint16_t cols;
	uint32_t buf_size;
	float *buf;
}mat_f32s;

typedef struct
{
	enum mat_type type;
	uint16_t rows;
	uint16_t cols;
	uint32_t buf_size;
	int16_t *buf;
}mat_s16_s;

mat_s *mat_clone_u8_f32(mat_s *src);

mat_s *mat_create(uint16_t rows, uint16_t cols, enum mat_type type);
void mat_free(mat_s *mat);
void mat_set_zero_f32(mat_s *mat);
void mat_set_zero_f64(mat_s *mat);
void mat_convert_u8_to_f32(mat_s *src, mat_s *dst);
mat_s *mat_resize_half_f32(mat_s *src);
void mat_resize_double_f32(mat_s *src, mat_s *dst);
mat_s *mat_convert_f32_scale(mat_s *src, float scale);
//mat_s *mat_convert_f32_to_u8(mat_s *src);
void mat_convert_f32_to_u8(mat_s *src, mat_s *dst);
void mat_printf_u8(uint8_t *mat, int m, int n);
void mat_printf_u16(uint16_t *mat, int m, int n);
void mat_printf_u32(uint32_t *mat, int m, int n);
void mat_printf_f32(float *mat, int m, int n);
void mat_printf_f64(double *mat, int m, int n);
mat_s *mat_clone_f32(mat_s *src);
mat_s *mat_from_data_u8(uint8_t *data, int rows, int cols);
//mat_s *mat_convert_u8_to_f32(mat_s *src);

void write_mat_data_to_file_u8(char *str, mat_s *mat);
void write_mat_data_to_file_f32(char *str, mat_s *mat);

void mat_set_f32(mat_s* mat, int r, int c, float val);
void mat_set_f64(mat_s* mat, int r, int c, double val);
mat_s *mat_clone_f32_test(mat_s *src);
mat_s *mat_clone_u8_test(mat_s *src);


/**/
mat_s* u8_data_to_mat_u8(uint8_t* data, int rows, int cols);
mat_s* array_f32_to_mat_f32(float *array, int rows, int cols);
mat_s* array_f64_to_mat_f32(double *array, int rows, int cols);
mat_s* array_f64_to_mat_f64(double *array, int rows, int cols);
double* mat_f32_to_array_f64(mat_s* mat);
double* mat_f64_to_array_f64(mat_s* mat);
float* mat_f32_to_array_f32(mat_s* mat);
mat_s* mat_get_row(mat_s* mat,int row);
mat_s* mat_get_col(mat_s* mat,int col);
double *mat_get_row(double* mat, int rows,int cols,int row);
double *mat_get_col(double* mat, int rows,int cols,int col);

/**/

#endif
