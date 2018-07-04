#include "file_operation.h"
#include "convert.h"
#include <opencv2/highgui/highgui.hpp>

CvMat *mat_to_CvMat(mat_s *mat)
{
	CvMat *Mat = NULL;
	if (mat->type == MAT_U8)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_8UC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_8UC1, (uint8_t *)mat->buf);
	}
	else if (mat->type == MAT_S8)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_8SC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_8SC1, (int8_t *)mat->buf);
	}
	else if (mat->type == MAT_U16)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_16UC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_16UC1, (uint16_t *)mat->buf);
	}
	else if (mat->type == MAT_S16)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_16SC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_16SC1, (int16_t *)mat->buf);
	}
	else if (mat->type == MAT_F32)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_32FC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_32FC1, (float *)mat->buf);
	}
	else if (mat->type == MAT_F64)
	{
		Mat = cvCreateMat(mat->rows, mat->cols, CV_64FC1);
		cvInitMatHeader(Mat, mat->rows, mat->cols, CV_64FC1, (double *)mat->buf);
	}
	return Mat;
}

CvMat data_to_CvMat(uint8_t *data, int rows, int cols)
{
	CvMat Mat = cvMat(rows,cols,CV_8UC1,data);
	return Mat;
}

IplImage *mat_to_image(mat_s *mat)
{
	CvMat Mat;
	if (mat->type = MAT_U8)
	{
		Mat = cvMat(mat->rows, mat->cols, CV_8UC1, mat->buf);
	}
	else if (mat->type = MAT_S8)
	{
		Mat = cvMat(mat->rows, mat->cols, CV_8SC1, mat->buf);
	}
	else if (mat->type = MAT_U16)
	{
		Mat = cvMat(mat->rows, mat->cols, CV_16UC1, mat->buf);
	}
	else if (mat->type = MAT_S16)
	{
		Mat = cvMat(mat->rows, mat->cols, CV_16SC1, mat->buf);
	}
	else if (mat->type = MAT_F32)
	{
		Mat = cvMat(mat->rows, mat->cols, CV_32FC1, mat->buf);
	}
	
	IplImage *img = cvCreateImage(cvGetSize(&Mat), 8, 1);
	cvGetImage(&Mat, img);
	return img;
}

IplImage *data_to_image(uint8_t *data, int rows, int cols)
{
	CvMat Mat = data_to_CvMat(data, rows, cols);
	IplImage *img = cvCreateImage(cvGetSize(&Mat), 8, 1);
	cvGetImage(&Mat, img);
	return img;
}

void image_show_mat(char *str,mat_s *mat)
{
	IplImage *img = mat_to_image(mat);
	cvNamedWindow(str, 1);
	cvShowImage(str, img);
	cvWaitKey(0);
	cvDestroyWindow(str);
//	cvReleaseImage(&img);
}

void image_show_data(char *str,uint8_t *data, int rows, int cols)
{
	IplImage *img = data_to_image(data,rows,cols);
	cvNamedWindow(str, 1);
	cvShowImage(str, img);
	cvWaitKey(0);
	cvDestroyWindow(str);
//	cvReleaseImage(&img);
}


void test_convert()
{
	char s[] = "pic1";
	IplImage *img = NULL;
	uint8_t *p;
	p = (uint8_t *)malloc(480 * 640);
	if (p) 
	{
	}
	else 
	{
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
	image_show_data(s,p, m, n);
	free(p);
}



