#ifndef _CONVERT_H_
#define _CONVERT_H_
#include <stdio.h>
#include <stdint.h>
#include <opencv\cv.h>
#include "mat.h" 

CvMat *mat_to_CvMat(mat_s *mat);
CvMat data_to_CvMat(uint8_t *data, int rows, int cols);
IplImage *mat_to_image(mat_s *mat);
IplImage *data_to_image(uint8_t *data, int rows, int cols);
extern "C" void image_show_mat(char *str, mat_s *mat);
extern "C" void image_show_data(char *str, uint8_t *mat, int rows, int cols);
void test_convert();

#endif
