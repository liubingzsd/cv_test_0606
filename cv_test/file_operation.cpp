#include <malloc.h>
#include "file_operation.h"

void read_data_from_file_f32(FILE *fp, float *data, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf_s(fp, "%f ", data + i * n + j);
		}
	}
}
void read_data_from_file_u8(FILE *fp, uint8_t *data, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fscanf_s(fp, "%hhd", data + i * n + j);
		}
	}
}

void write_data_to_file_f32(FILE *fp, float *data, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fprintf(fp, "%f ", *(data + i * n + j));
		}
		fprintf(fp, "\n");
	}
}

void write_data_to_file_u8(FILE *fp, uint8_t *data, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fprintf(fp, "%d ", *(data + i * n + j));
		}
		fprintf(fp, "\n");
	}
}

void write_data_to_file_s16(FILE *fp, int16_t *data, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			fprintf(fp, "%d ", *(data + i * n + j));
		}
		fprintf(fp, "\n");
	}
}

mat_s* mat_read_data_u8_from_file_path(char* path,int rows,int cols)
{
	FILE* fp = fopen(path, "r");
	mat_s* mat = NULL;
	if (fp == NULL)
	{
		printf("NO file fp\n");
	}
	else
	{
		mat = mat_create(rows, cols, MAT_U8);
		read_data_from_file_u8(fp, (uint8_t*)mat->buf, rows, cols);
	}
	return mat;
}

uint8_t* read_data_u8_from_file_path(char* path, int rows, int cols)
{
	FILE* fp = fopen(path, "r");
	uint8_t* data = NULL;
	if (fp == NULL)
	{
		printf("NO file fp\n");
	}
	else
	{
		data = (uint8_t *)malloc(sizeof(uint8_t) * rows * cols);
		read_data_from_file_u8(fp, data, rows, cols);
	}
	return data;
}


void test_file_wr()
{
	uint8_t data[480 * 640] = { 0 };
	uint8_t data1[480 * 640] = { 0 };
	uint8_t data2[640] = { 0 };
	int m = 480;
	int n = 640;
	FILE *fp = fopen("data.dat", "r");
	if (fp == NULL)
	{
		printf("NO file fp\n");
	}
	else
	{
		read_data_from_file_u8(fp, data, m, n);
	}

	FILE *fp1 = fopen("data1.dat", "w");
	if (fp == NULL)
	{
		printf("NO file fp1\n");
	}
	else
	{
		write_data_to_file_u8(fp1, data, m, n);
	}
	fclose(fp1);
}


void test_read_file_u8()
{
	uint8_t *p;
	p = (uint8_t *)malloc(480 * 640);
	if (p) {
	}
	else {
		printf("mem allocated filed \n");
	}
	int m = 480;
	int n = 640;
	FILE *fp = fopen("data.dat", "r");
	if (fp == NULL)
	{
		printf("NO file fp\n");
	}
	else
	{
		read_data_from_file_u8(fp, p, m, n);
	}
	fclose(fp);
	FILE *fp1 = fopen("data1.dat", "w+");
	if (fp1 == NULL)
	{
		printf("NO file fp1\n");
	}
	else
	{
		write_data_to_file_u8(fp1, p, m, n);
	}
	fclose(fp1);
	free(p);
}
void test_file_read_float()
{
	float data[400][2] = { 0 };
	int m = 400;
	int n = 2;
	FILE *fp = fopen("data.dat", "r");
	if (fp == NULL)
	{
		printf("NO file\n");
	}
	else
	{
		read_data_from_file_f32(fp, (float *)data, m, n);
	}
}
