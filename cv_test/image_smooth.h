#ifndef _IMAGE_SMOOTH_H_
#define _IMAGE_SMOOTH_H_

void image_smooth_gaussian(mat_s *src, mat_s *dst,float sig1, float sig2);
mat_s* image_smooth_gaussian_rt(mat_s *src, float sig1, float sig2);
#endif
