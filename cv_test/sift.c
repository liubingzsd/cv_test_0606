#include "mat.h"
#include "sift.h"
#include "kernel.h"
#include "common.h"
#include "mat_math.h"
#include "file_operation.h"
#include "image_filter.h"

#define GAUSSKERN				3.5f
#define SCALESPEROCTAVE			2
#define MAXOCTAVES				4
#define CURVATURE_THRESHOLD     10.0f
#define DOUBLE_BASE_IMAGE_SIZE  1
#define CONTRAST_THRESHOLD      0.02f 
#define PEAK_REAL_THRESH        0.8f
#define GRID_SPACING            4 

image_octave_s *dog_octaves, *mag_pyr, *grad_pyr;
key_point_s *key_points = NULL;
key_point_s *key_descriptors = NULL;

void image_show_data(char *str, uint8_t *data, int rows, int cols);
void image_show_mat(char *str, mat_s *mat);

void __image_show_f32(char *str1, mat_s *mat, int one_sign)
{
	mat_s *mat_scale = NULL, *mat_u8 = NULL;
	if (one_sign)
	{
		mat_scale = mat_create(mat->rows, mat->cols, MAT_F32);
		mat_u8 = mat_create(mat->rows, mat->cols, MAT_U8);
		mat_scale = mat_convert_f32_scale(mat, 255.0f);
		mat_convert_f32_to_u8(mat_scale, mat_u8);
		write_mat_data_to_file_u8("mat_u8",mat_u8);
		image_show_mat(str1, mat_u8);
	}
	else
	{
		mat_u8 = mat_create(mat->rows, mat->cols, MAT_U8);
		mat_convert_f32_to_u8(mat, mat_u8);
		image_show_mat(str1, mat_u8);
	}
}

void image_show_u8(char *str1, mat_s *mat)
{
	image_show_mat(str1, mat);
}

uint8_t image_blur_ker_size(float sigma)
{
	uint8_t size;
	size = size = (uint8_t)max(3.0f, 2.0f * GAUSSKERN * sigma + 1.0f);
	if (size % 2 == 0)
	{
		size += 1;
	}
	return size;
}
 

void image_blur_f32(mat_s *src, mat_s *dst,float sigma)
{
	kernel_s ker;
	uint8_t size;
	mat_s *tmp;
	tmp = mat_create(src->rows, src->cols, MAT_F32);
	size = image_blur_ker_size(sigma);
	kernel_create_1d(&ker, size, KERNEL_F32);
	gauss_kernel_1d(&ker, sigma);
	mat_clone_f32(src,tmp);
	image_filter_1d_x(src, tmp, &ker, NONE);
	mat_clone_f32(tmp,dst);
	image_filter_1d_y(tmp, dst, &ker, NONE);
	mat_free(tmp);
}


/*
mat_s *image_blur_f32_rt(mat_s *src, float sigma)
{
	kernel_s ker;
	uint8_t size;
	mat_s *tmp = NULL,*dst = NULL;
	size = image_blur_ker_size(sigma);
	kernel_create_1d(&ker, size, KERNEL_F32);
	gauss_kernel_1d(&ker, sigma);
	tmp = mat_clone_f32_test(src);
	image_filter_1d_x(src, tmp, &ker, NONE);
	write_mat_data_to_file_f32("tmp.dat", tmp);
	dst = mat_clone_f32_test(tmp);
	write_mat_data_to_file_f32("dst.dat", dst);
	image_filter_1d_y(tmp, dst, &ker, NONE);
	write_mat_data_to_file_f32("dst_f.dat", dst);
//	__image_show_f32("dst", dst, 1);
	mat_free(tmp);
	return dst;
}
*/
mat_s *image_scale_init(mat_s *src)
{
	float init_sigma = 0.5f;
	float preblur_sigma = 1.0f;
	float sigma = 0.0f;
	mat_s *img_blur = mat_create(src->rows,src->cols,MAT_F32);
	mat_s *dst = NULL;
	image_blur_f32(src,img_blur,init_sigma);
	if (DOUBLE_BASE_IMAGE_SIZE)
	{
		uint32_t rows, cols;
		rows = 2 * src->rows;
		cols = 2 * src->cols;
		mat_s *img_double, *img_double_blur;
		img_double = mat_create(rows, cols, MAT_F32);
		dst = mat_create(rows, cols, MAT_F32);
		img_double_blur = mat_create(rows, cols, MAT_F32);
		mat_resize_double_f32(img_blur, img_double);
		image_blur_f32(img_double, img_double_blur,preblur_sigma);
		sigma = (float)sqrt((4.0f * init_sigma * init_sigma) + preblur_sigma * preblur_sigma);
		image_blur_f32(img_double_blur,dst,sigma);
		mat_free(img_blur);
		mat_free(img_double);
		mat_free(img_double_blur);
	}
	else
	{
		dst = mat_create(src->rows,src->cols,MAT_F32);
		sigma = (float)sqrt((4.0f * init_sigma * init_sigma) + preblur_sigma * preblur_sigma);
		image_blur_f32(img_blur,dst,sigma);
		mat_free(img_blur);
		__image_show_f32("dst", dst, 1);
	}
	return dst;
}

/*
image_octave_s *build_gaussian_octaves(mat_s *src)
{
	image_octave_s *octaves;
	mat_s *img_tmp, *dst, *temp;
	int i, j;
	float k;
	float initial_sigma, sigma, absolue_sigma, sigma_f;
	k = (float)pow(2.0f, 1.0f / (float)(SCALESPEROCTAVE));
	int dim = min(src->rows, src->cols);
	int num_octaves = (int)(log(dim) / log(2.0)) - 2;
	num_octaves = min(num_octaves, MAXOCTAVES);
	octaves = (image_octave_s *)malloc(num_octaves * sizeof(image_octave_s));
	dog_octaves = (image_octave_s *)malloc(num_octaves * sizeof(image_octave_s));
	img_tmp = mat_clone_f32(src);
	initial_sigma = (float)sqrt(2.0f);
	for (i = 0; i < num_octaves; i++)
	{
		printf("Building octave %d of dimesion (%d, %d) \n", i, img_tmp->rows, img_tmp->cols);
		octaves[i].octave = (image_level_s *)malloc((SCALESPEROCTAVE + 3) * sizeof(image_level_s));
		dog_octaves[i].octave = (image_level_s *)malloc((SCALESPEROCTAVE + 2) * sizeof(image_level_s));
		(octaves[i].octave)[0].level = img_tmp;
		octaves[i].cols = img_tmp->cols;
		octaves[i].rows = img_tmp->rows;
		dog_octaves[i].cols = img_tmp->cols;
		dog_octaves[i].rows = img_tmp->rows;
		if (DOUBLE_BASE_IMAGE_SIZE)
		{
			octaves[i].sub_sample = (float)pow(2.0f, i)*0.5f;
		}
		else
		{
			octaves[i].sub_sample = (float)pow(2.0f, i);
		}
		if (i == 0)
		{
			(octaves[0].octave)[0].level_sigma = initial_sigma;
			(octaves[0].octave)[0].absolute_sigma = initial_sigma;
			printf("%d scale and blur sigma : %f \n", i, (octaves[0].sub_sample) *((octaves[0].octave)[0].absolute_sigma));
		}
		else
		{
			(octaves[i].octave)[0].level_sigma = (octaves[i - 1].octave)[SCALESPEROCTAVE].level_sigma;
			(octaves[i].octave)[0].absolute_sigma = (octaves[i - 1].octave)[SCALESPEROCTAVE].absolute_sigma;
			printf("%d scale and blur sigma : %f \n", i, ((octaves[i].octave)[0].absolute_sigma));

		}
		sigma = initial_sigma;
		for (j = 1; j < SCALESPEROCTAVE + 3; j++)
		{
			dst = mat_create(img_tmp->rows, img_tmp->cols, MAT_F32);
			temp = mat_create(img_tmp->rows, img_tmp->cols, MAT_F32);
			sigma_f = (float)sqrt(k * k - 1)*sigma;
			sigma = k * sigma;
			absolue_sigma = sigma * (octaves[i].sub_sample);
			printf("%d scale and blur_sigma %f \n", j, absolue_sigma);
			(octaves[i].octave)[j].level_sigma = sigma;
			(octaves[i].octave)[j].absolute_sigma = absolue_sigma;
			int length = image_blur_ker_size(sigma_f);
			image_blur_f32((octaves[i].octave)[j - 1].level, dst, sigma_f);
			(octaves[i].octave)[j].level_sigma_length = length;
			(octaves[i].octave)[j].level = dst;
			temp = mat_sub((octaves[i].octave)[j].level, (octaves[i].octave)[j - 1].level);
			(dog_octaves[i].octave)[j - 1].level = temp;
		}
		img_tmp = mat_resize_half_f32(((octaves[i].octave)[SCALESPEROCTAVE].level));
	}
	return octaves;
}

float img_level(uint8_t i, uint8_t j, uint16_t m, uint16_t n)
{
	float tmp;
	float *buf = (float *)(dog_octaves[i]).octave[j].level->buf;
	uint16_t col = (dog_octaves[i]).octave[j].level->cols;
	tmp = buf[m * col + n];
	return tmp;
}


int detect_key_point(int num_octaves, image_octave_s *gaussian_pyr)
{
	float curvature_threshold;
	curvature_threshold = (float)((CURVATURE_THRESHOLD + 1)*(CURVATURE_THRESHOLD + 1)) / CURVATURE_THRESHOLD;
	int key_point_count = 0;
	for (int i = 0; i < num_octaves; i++)
	{
		for (int j = 1; j < SCALESPEROCTAVE + 1; j++)
		{
			int dim = (int)(0.5f * (gaussian_pyr[i].octave)[j].level_sigma_length + 0.5f);
			for (int m = dim; m < dog_octaves[i].rows - dim; m++)
			{
				for (int n = dim; n < dog_octaves[i].cols - dim; n++)
				{
					if (fabs(img_level(i, j, m, n)) >= CONTRAST_THRESHOLD)
					{
						if (img_level(i, j, m, n) != 0.0f)
						{
							float inf_val = img_level(i, j, m, n);
							if (((inf_val <= img_level(i, j - 1, m - 1, n - 1)) &&
								(inf_val <= img_level(i, j - 1, m, n - 1)) &&
								(inf_val <= img_level(i, j - 1, m + 1, n - 1)) &&
								(inf_val <= img_level(i, j - 1, m - 1, n)) &&
								(inf_val <= img_level(i, j - 1, m, n)) &&
								(inf_val <= img_level(i, j - 1, m + 1, n)) &&
								(inf_val <= img_level(i, j - 1, m - 1, n + 1)) &&
								(inf_val <= img_level(i, j - 1, m, n + 1)) &&
								(inf_val <= img_level(i, j - 1, m + 1, n + 1)) &&    //底层的小尺度9

								(inf_val <= img_level(i, j, m - 1, n - 1)) &&
								(inf_val <= img_level(i, j, m, n - 1)) &&
								(inf_val <= img_level(i, j, m + 1, n - 1)) &&
								(inf_val <= img_level(i, j, m - 1, n)) &&
								(inf_val <= img_level(i, j, m + 1, n)) &&
								(inf_val <= img_level(i, j, m - 1, n + 1)) &&
								(inf_val <= img_level(i, j, m, n + 1)) &&
								(inf_val <= img_level(i, j, m + 1, n + 1)) &&     //当前层8

								(inf_val <= img_level(i, j + 1, m - 1, n - 1)) &&
								(inf_val <= img_level(i, j + 1, m, n - 1)) &&
								(inf_val <= img_level(i, j + 1, m + 1, n - 1)) &&
								(inf_val <= img_level(i, j + 1, m - 1, n)) &&
								(inf_val <= img_level(i, j + 1, m, n)) &&
								(inf_val <= img_level(i, j + 1, m + 1, n)) &&
								(inf_val <= img_level(i, j + 1, m - 1, n + 1)) &&
								(inf_val <= img_level(i, j + 1, m, n + 1)) &&
								(inf_val <= img_level(i, j + 1, m + 1, n + 1))     //下一层大尺度9        
								) ||
								((inf_val >= img_level(i, j - 1, m - 1, n - 1)) &&
								(inf_val >= img_level(i, j - 1, m, n - 1)) &&
									(inf_val >= img_level(i, j - 1, m + 1, n - 1)) &&
									(inf_val >= img_level(i, j - 1, m - 1, n)) &&
									(inf_val >= img_level(i, j - 1, m, n)) &&
									(inf_val >= img_level(i, j - 1, m + 1, n)) &&
									(inf_val >= img_level(i, j - 1, m - 1, n + 1)) &&
									(inf_val >= img_level(i, j - 1, m, n + 1)) &&
									(inf_val >= img_level(i, j - 1, m + 1, n + 1)) &&

									(inf_val >= img_level(i, j, m - 1, n - 1)) &&
									(inf_val >= img_level(i, j, m, n - 1)) &&
									(inf_val >= img_level(i, j, m + 1, n - 1)) &&
									(inf_val >= img_level(i, j, m - 1, n)) &&
									(inf_val >= img_level(i, j, m + 1, n)) &&
									(inf_val >= img_level(i, j, m - 1, n + 1)) &&
									(inf_val >= img_level(i, j, m, n + 1)) &&
									(inf_val >= img_level(i, j, m + 1, n + 1)) &&

									(inf_val >= img_level(i, j + 1, m - 1, n - 1)) &&
									(inf_val >= img_level(i, j + 1, m, n - 1)) &&
									(inf_val >= img_level(i, j + 1, m + 1, n - 1)) &&
									(inf_val >= img_level(i, j + 1, m - 1, n)) &&
									(inf_val >= img_level(i, j + 1, m, n)) &&
									(inf_val >= img_level(i, j + 1, m + 1, n)) &&
									(inf_val >= img_level(i, j + 1, m - 1, n + 1)) &&
									(inf_val >= img_level(i, j + 1, m, n + 1)) &&
									(inf_val >= img_level(i, j + 1, m + 1, n + 1))
									))      //2、满足26个中极值点
							{
								if (fabs(img_level(i, j, m, n)) >= CONTRAST_THRESHOLD)
								{
									float dxx, dyy, dxy, tr_h, det_h, curvature_ratio;
									dxx = img_level(i, j, m, n - 1) + img_level(i, j, m, n + 1) - 2.0f * img_level(i, j, m, n);
									dyy = img_level(i, j, m - 1, n) + img_level(i, j, m + 1, n) - 2.0f * img_level(i, j, m, n);
									dxy = img_level(i, j, m - 1, n - 1) + img_level(i, j, m + 1, n + 1) - img_level(i, j, m + 1, n - 1) - img_level(i, j, m - 1, n + 1);
									tr_h = dxx + dyy;
									det_h = dxx * dyy - dxy * dxy;
									curvature_ratio = (1.0f * tr_h * tr_h) / det_h;
									if (det_h > 0.0f && curvature_ratio <= curvature_threshold)
									{
										key_point_count++;
										key_point_s *k;
										k = (key_point_s *)malloc(sizeof(key_point_s));
										k->next = key_points;
										key_points = k;
										k->row = m * gaussian_pyr[i].sub_sample;
										k->col = n * gaussian_pyr[i].sub_sample;
										k->sy = (float)m;
										k->sx = (float)n;
										k->octave = i;
										k->level = j;
										k->scale = (gaussian_pyr[i].octave)[j].absolute_sigma;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return key_point_count;
}

float mat_val(mat_s *mat, uint16_t i, uint16_t j)
{
	float *mat_buf = (float *)mat->buf;
	return mat_buf[i * mat->cols + j];
}


void calc_grad_direc_mag(int num_octaves, image_octave_s *gaussian_pyr)
{
	float mat1, mat2, ori_val;
	mag_pyr = (image_octave_s *)malloc(num_octaves * sizeof(image_octave_s));
	grad_pyr = (image_octave_s *)malloc(num_octaves * sizeof(image_octave_s));
	for (int i = 0; i < num_octaves; i++)
	{
		mag_pyr[i].octave = (image_level_s *)malloc(SCALESPEROCTAVE * sizeof(image_level_s));
		grad_pyr[i].octave = (image_level_s *)malloc(SCALESPEROCTAVE * sizeof(image_level_s));
		for (int j = 1; j < SCALESPEROCTAVE + 1; j++)
		{
			mat_s *mag = mat_create(gaussian_pyr[i].rows, gaussian_pyr[i].cols, MAT_F32);
			mat_s *ori = mat_create(gaussian_pyr[i].rows, gaussian_pyr[i].cols, MAT_F32);
			mat_s *temp_mat1 = mat_create(gaussian_pyr[i].rows, gaussian_pyr[i].cols, MAT_F32);
			mat_s *temp_mat2 = mat_create(gaussian_pyr[i].rows, gaussian_pyr[i].cols, MAT_F32);
			mat_set_zero_f32(mag);
			mat_set_zero_f32(ori);
			mat_set_zero_f32(temp_mat1);
			mat_set_zero_f32(temp_mat2);
			float *mat1_buf = (float *)temp_mat1->buf;
			float *mat2_buf = (float *)temp_mat2->buf;
			float *mag_buf = (float *)mag->buf;
			float *ori_buf = (float *)ori->buf;
			for (int m = 1; m < (gaussian_pyr[i].rows - 1); m++)
			{
				for (int n = 1; n < (gaussian_pyr[i].cols - 1); n++)
				{
					mat1 = mat1_buf[m * temp_mat1->cols + n] = 0.5f * (img_level(i, j, m, n + 1) - img_level(i, j, m, n - 1));
					mat2 = mat2_buf[m * temp_mat2->cols + n] = 0.5f * (img_level(i, j, m + 1, n) - img_level(i, j, m - 1, n));
					mag_buf[m * mag->cols + n] = (float)sqrt(mat1 * mat1 + mat2 * mat2);
					ori_val = ori_buf[m * ori->cols + n] = (float)atan(mat1 / mat2);
					if (ori_val == PI)
					{
						ori_val = -PI;
					}
				}
			}
			(mag_pyr[i].octave)[j - 1].level = mag;
			(grad_pyr[i].octave)[j - 1].level = ori;
			mat_free(temp_mat1);
			mat_free(temp_mat2);
		}
	}
}

int find_closest_rotation_bin(int bin_count, float angle)
{
	angle += PI;
	angle /= 2.0f * PI;
	angle *= bin_count;
	int idx = (int)angle;
	if (idx == bin_count)
	{
		idx = 0;
	}
	return idx;
}

void average_weak_bins(float *hist, int bin_count)
{
	for (int sn = 0; sn < 2; sn++)
	{
		float first_e = hist[0];
		float last = hist[bin_count - 1];
		for (int sw = 0; sw < bin_count; sw++)
		{
			float cur = hist[sw];
			float next = (sw == (bin_count - 1)) ? first_e : hist[(sw + 1) % bin_count];
			hist[sw] = (last + cur + next) / 3.0f;
			last = cur;
		}
	}
}

bool inter_polate_orientation(float left, float middle, float right, float *degree_correction, float *peak_value)
{
	float a = ((left + right) - 2.0f * middle) / 2.0f;
	if (a == 0.0f)
	{
		return false;
	}
	float c = (((left - middle) / a) - 1.0f) / 2.0f;
	float b = middle - c * c * a;
	if (c < -0.5f || c > 0.5f)
	{
		return false;
	}
	*degree_correction = c;
	*peak_value = b;
	return true;
}

void assign_the_main_orientation(int num_octaves, image_octave_s *gaussian_pyr, image_octave_s *mag_pyr, image_octave_s *grad_pyr)
{
	int num_bins = 36;
	float hist_step = 2.0f * PI / num_bins;
	float hist_orient[36] = { 0.0f };
	for (int i = 0; i < 36; i++)
	{
		hist_orient[i] = -PI + i * hist_step;
	}
	float sigma1 = ((gaussian_pyr[0].octave)[SCALESPEROCTAVE].absolute_sigma) / (gaussian_pyr[0].sub_sample);
	int zero_pad = (int)(max(3.0f, 2 * GAUSSKERN*sigma1 + 1.0f)*0.5f + 0.5f);
	int key_point_count = 0;
	key_point_s *p = key_points;
	while (p)
	{
		int i = p->octave;
		int j = p->level;
		int m = (int)p->sy;
		int n = (int)p->sx;
		if ((m > zero_pad) && (m < gaussian_pyr[i].rows - zero_pad) && (n >= zero_pad) && (n < gaussian_pyr[i].cols - zero_pad))
		{
			float sigma = ((gaussian_pyr[i].octave)[j].absolute_sigma) / (gaussian_pyr[i].sub_sample);
			mat_s *mat = gaussian_kernel_2d(sigma);
			float *mat_buf = (float *)mat->buf;
			int dim = (int)(0.5f * (mat->rows));
			float *orient_hist = (float *)malloc(36 * sizeof(float));
			for (int sw = 0; sw < 36; sw++)
			{
				orient_hist[sw] = 0;
			}
			for (int x = m - dim, mm = 0; x <= (m + dim); x++, mm++)
			{
				for (int y = n - dim, nn = 0; y <= (n + dim); y++, nn++)
				{
					float dx = 0.5f * img_level(i, j, x, y + 1) - img_level(i, j, x, y - 1);
					float dy = 0.5f * img_level(i, j, x + 1, y) - img_level(i, j, x - 1, y);
					float mag = (float)sqrt(dx * dx + dy * dy);
					float ori = (float)atan(1.0f * dy / dx);
					int bin_idx = find_closest_rotation_bin(36, ori);
					orient_hist[bin_idx] = orient_hist[bin_idx] + 1.0f * mag * mat_buf[mm * mat->cols + nn];
				}
			}
			average_weak_bins(orient_hist, 36);
			float max_grad = 0.0f;
			int max_bin = 0;
			for (int b = 0; b < 36; b++)
			{
				if (orient_hist[b] > max_grad)
				{
					max_grad = orient_hist[b];
					max_bin = b;
				}
			}
			float max_peak_value = 0.0f;
			float max_degree_correction = 0.0f;
			if ((inter_polate_orientation(orient_hist[max_bin] == 0 ? (36 - 1) : (max_bin - 1),orient_hist[max_bin], orient_hist[(max_bin + 1) % 36],
				&max_degree_correction, &max_peak_value)) == false)
			{
				printf("bug parabola filtting broken \n");
			}
			bool bin_is_key_point[36];
			for (int b = 0; b < 36; b++)
			{
				bin_is_key_point[b] = false;
				if (b == max_bin)
				{
					bin_is_key_point[b] = true;
					continue;
				}
				if (orient_hist[b] < (PEAK_REAL_THRESH * max_peak_value))
				{
					continue;
				}
				int left_i = (b == 0) ? (36 - 1) : (b - 1);
				int right_i = (b + 1) % 36;
				if (orient_hist[b] <= orient_hist[left_i] || orient_hist[b] <= orient_hist[right_i])
				{
					continue;
				}
				bin_is_key_point[b] = true;
			}

			float one_bin_rad = (2.0f * PI) / 36;
			for (int b = 0; b < 36; b++)
			{
				if (bin_is_key_point[b] == false)
				{
					continue;
				}
				int b_left = (b == 0) ? (36 - 1) : (b - 1);
				int b_right = (b + 1) % 36;
				float peak_value, degree_correction;
				if (inter_polate_orientation(orient_hist[max_bin == 0 ? (36 - 1) : (max_bin - 1)],
					orient_hist[max_bin], orient_hist[(max_bin + 1) % 36],
					&degree_correction, &peak_value) == false)
				{
					printf("bug : parabola fitting broken \n");
				}
				float degree = (b + degree_correction) * one_bin_rad - PI;
				if (degree < -PI)
				{
					degree += 2.0f * PI;
				}
				else if (degree > PI)
				{
					degree -= 2.0f * PI;
				}


				key_point_s *k;
				k = (key_point_s *)malloc(sizeof(key_point_s));
				k->next = key_descriptors;
				key_descriptors = k;
				k->descriptor = (float *)malloc(128 * sizeof(float));
				k->row = p->row;
				k->sy = p->sy;
				k->sx = p->sx;
				k->octave = p->octave;
				k->level = p->level;
				k->scale = (float)p->level;
				k->scale = p->scale;
				k->ori = degree;
				k->mag = peak_value;
			}
			free(orient_hist);
		}
		p = p->next;
	}
}


float get_pixel_bilinear(mat_s *im, float col, float row)
{
	int irow, icol;
	float rfrac, cfrac;
	float row1 = 0, row2 = 0;
	int width = im->cols;
	int height = im->rows;
	irow = (int)row;
	icol = (int)col;
	if (irow < 0 || irow >= height || icol < 0 || icol >= width)
	{
		return 0;
	}
	if (row > height - 1)
	{
		row = (float)(height - 1);
	}
	if (col > width - 1)
	{
		col = (float)(width - 1);
	}
	rfrac = 1.0f - (row - (float)irow);
	cfrac = 1.0f - (col - (float)icol);
	if (cfrac < 1)
	{
		row1 = cfrac * mat_val(im, irow, icol) + (1.0f - cfrac) * mat_val(im, irow, icol + 1);
	}
	else
	{
		row1 = mat_val(im, irow, icol);
	}
	if (rfrac < 1)
	{
		if (cfrac < 1)
		{
			row2 = cfrac * mat_val(im, irow + 1, icol) + (1.0f - cfrac) * mat_val(im, irow + 1, icol + 1);
		}
		else
		{
			row2 = mat_val(im, irow + 1, icol);
		}
	}
	return rfrac * row1 + (1.0f - rfrac) * row2;
}

float get_vec_norm(float *vec, int dim)
{
	float sum = 0.0f;
	for (int i = 0; i < dim; i++)
	{
		sum += vec[i] * vec[i];
	}
	return (float)sqrt(sum);
}

void extract_feature_descriptors(int num_octaves, image_octave_s *gaussian_pyr)
{
	float orient_bin_spacing = PI / 4.0f;
	float orient_angles[8] = { -PI,-PI + orient_bin_spacing,-PI * 0.5f,-orient_bin_spacing,0.0f,orient_bin_spacing,PI*0.5f,PI + orient_bin_spacing };
	float *feat_grid = (float *)malloc(2 * 16 * sizeof(float));
	for (int i = 0; i < GRID_SPACING; i++)
	{
		for (int j = 0; j < 2 * GRID_SPACING; ++j, ++j)
		{
			feat_grid[i * 2 * GRID_SPACING + j] = -6.0f + i * GRID_SPACING;
			feat_grid[i * 2 * GRID_SPACING + j + 1] = -6.0f + 0.5f * j * GRID_SPACING;
		}
	}

	float *feat_samples = (float *)malloc(2 * 256 * sizeof(float));
	for (int i = 0; i < 4 * GRID_SPACING; i++)
	{
		for (int j = 0; j < 8 * GRID_SPACING; j += 2)
		{
			feat_samples[i * 8 * GRID_SPACING + j] = -(2 * GRID_SPACING - 0.5f) + i;
			feat_samples[i * 8 * GRID_SPACING + j + 1] = -(2 * GRID_SPACING - 0.5f) + 0.5f * j;
		}
	}
	float feat_window = 2 * GRID_SPACING;
	key_point_s *p = key_descriptors;
	while (p)
	{
		float scale = (gaussian_pyr[p->octave].octave)[p->level].absolute_sigma;
		float sine = (float)sin(p->ori);
		float cosine = (float)cos(p->ori);
		float *feat_center = (float *)malloc(2 * 16 * sizeof(float));
		for (int i = 0; i < GRID_SPACING; i++)
		{
			for (int j = 0; j < 2 * GRID_SPACING; j += 2)
			{
				float x = feat_grid[i * 2 * GRID_SPACING + j];
				float y = feat_grid[i * 2 * GRID_SPACING + j + 1];
				feat_center[i * 2 * GRID_SPACING + j] = ((cosine * x + sine * y) + p->sx);
				feat_center[i * 2 * GRID_SPACING + j + 1] = ((-sine * x + cosine * y) + p->sy);
			}
			float *feat = (float *)malloc(2 * 256 * sizeof(float));
			for (int i = 0; i < 64 * GRID_SPACING; i++, i++)
			{
				float x = feat_samples[i];
				float y = feat_samples[i + 1];
				feat[i] = ((cosine * x + sine * y) + p->sx);
				feat[i + 1] = (-sine * x + cosine * y) + p->sy;
			}
			float * feat_desc = (float *)malloc(128 * sizeof(float));
			for (int i = 0; i < 128; i++)
			{
				feat_desc[i] = 0.0f;
			}
			for (int i = 0; i < 512; ++i, ++i)
			{
				float x_sample = feat[i];
				float y_sample = feat[i + 1];
				float sample12 = get_pixel_bilinear(((gaussian_pyr[p->octave].octave)[p->level]).level, x_sample, y_sample - 1);
				float sample21 = get_pixel_bilinear(((gaussian_pyr[p->octave].octave)[p->level]).level, x_sample - 1, y_sample);
				float sample22 = get_pixel_bilinear(((gaussian_pyr[p->octave].octave)[p->level]).level, x_sample, y_sample);
				float sample23 = get_pixel_bilinear(((gaussian_pyr[p->octave].octave)[p->level]).level, x_sample + 1, y_sample);
				float sample32 = get_pixel_bilinear(((gaussian_pyr[p->octave].octave)[p->level]).level, x_sample, y_sample + 1);
				float diff_x = sample23 - sample12;
				float diff_y = sample32 - sample12;
				float mag_sample = (float)sqrt(diff_x * diff_x + diff_y * diff_y);
				float grad_sample = (float)(atan(diff_y / diff_x));
				if (grad_sample == PI)
				{
					grad_sample = -PI;
				}
				float *x_weight = (float *)malloc(GRID_SPACING * GRID_SPACING * sizeof(float));
				float *y_weight = (float *)malloc(GRID_SPACING * GRID_SPACING * sizeof(float));
				float *pos_weight = (float *)malloc(8 * GRID_SPACING * GRID_SPACING * sizeof(float));
				for (int m = 0; m < 16; ++m)
				{
					for (int n = 0; n < 8; ++n)
					{
						pos_weight[m * 8 + n] = x_weight[m] * y_weight[m];
					}
				}
				free(x_weight);
				free(y_weight);
				float diff[8], orient_weight[128];
				for (int m = 0; m < 8; ++m)
				{
					float angle = (float)(grad_sample - (p->ori) - orient_angles[m] + PI);
					float temp = (float)angle / (2.0f * PI);
					angle -= (int)(temp) * (2.0f * PI);
					diff[m] = angle - PI;
				}
				float x = p->sx;
				float y = p->sy;
				float g = (float)exp(-((x_sample - x) * (x_sample - x) + (y_sample - y) * (y_sample - y)) / (2 * feat_window * feat_window)) / (2 * PI * feat_window);
				for (int m = 0; m < 128; ++m)
				{
					orient_weight[m] = (float)max((1.0f - 1.0f * fabs(diff[m % 8]) / orient_bin_spacing), 0);
					feat_desc[m] = feat_desc[m] + orient_weight[m] * pos_weight[m] * g * mag_sample;
				}
				free(pos_weight);
			}
			free(feat);
			free(feat_center);
			float norm = get_vec_norm(feat_desc, 128);
			for (int m = 0; m < 128; m++)
			{
				feat_desc[m] /= norm;
				if (feat_desc[m] > 0.2f)
				{
					feat_desc[m] = 0.2f;
				}
			}
			norm = get_vec_norm(feat_desc, 128);
			for (int m = 0; m < 128; m++)
			{
				feat_desc[m] /= norm;
				printf("%f ", feat_desc[m]);
			}
			printf("\n");
			p->descriptor = feat_desc;
			p = p->next;
		}
		free(feat_grid);
		free(feat_samples);
	}
}
*/

mat_s *read_image_data_u8()
{
	char str[] = "pic1";
	uint8_t *p;
	mat_s *mat = NULL;
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
	mat = mat_from_data_u8(p,m,n);
//	image_show_mat(str,mat);
	free(p);
	return mat;
}


void test_sift()
{
	int m = 480;
	int n = 640;
	int num_octaves = 0;
	image_octave_s *gaussian_pyr;
	mat_s *img_f32 = NULL,*img_u8 = NULL;
	mat_s *img_tmp = NULL;
	char str1[] = "test";
	mat_s *raw_img = NULL;
	mat_s *img_r_u8 = NULL, *img_r_f32 = NULL;
	raw_img = read_image_data_u8();
	img_r_u8 = mat_clone_u8_test(raw_img);
	img_f32 = mat_create(m,n,MAT_F32);
	mat_convert_u8_to_f32(raw_img,img_f32);
	img_r_f32 = mat_clone_f32_test(img_f32);

//	write_mat_data_to_file_f32("image_f32.dat",img_f32);
//	mat_convert_f32_to_u8(img_f32, img_u8);
//	write_mat_data_to_file_u8("image_u8.dat", img_u8);
//	__image_show_f32(str1, img_f32, 0);

	img_f32 = mat_convert_f32_scale(img_f32, 1.0f / 255);
	int dim = min(img_f32->rows, img_f32->cols);
	num_octaves = (int)(log((float)dim) / log(2.0f)) - 2;
	num_octaves = min(num_octaves, MAXOCTAVES);
	img_tmp = image_scale_init(img_f32);
	__image_show_f32("scale_img", img_tmp, 1);

/*
	gaussian_pyr = build_gaussian_octaves(img_tmp);
	int key_point_num = detect_key_point(num_octaves, gaussian_pyr);
	calc_grad_direc_mag(num_octaves, gaussian_pyr);
	assign_the_main_orientation(num_octaves, gaussian_pyr, mag_pyr, grad_pyr);
	extract_feature_descriptors(num_octaves, gaussian_pyr);
*/
}