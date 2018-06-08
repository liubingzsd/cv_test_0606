#ifndef _SIFT_TEST_H_
#define _SIFT_TEST_H_
#include <stdbool.h>
#include "mat.h"

typedef struct
{
	float level_sigma;
	float absolute_sigma;
	int level_sigma_length;
	mat_s *level;
}image_level_s;

typedef struct
{
	int rows, cols;
	float sub_sample;
	image_level_s *octave;
}image_octave_s;

typedef struct key_point
{
	float row, col;  /* keypoints in raw image's location*/
	float sx, sy;    /* keypoints in image_pyr's location*/
	int octave, level;/*keypoints octave and level*/
	float scale, ori, mag; /**/
	float *descriptor;    /**/
	struct key_point *next;
}key_point_s;


void test_sift();

#endif