#include <math.h>
#include <string.h>
#include "common.h"
#include "mat.h"
#include "svd_solver.h"
#include "convert.h"
#include "mat_math.h"
#include "sift_test.h"
#include "image_smooth.h"
#include "image_features.h"
#include "file_operation.h"
#include "feature_seq.h"
#include "kdtree.h"
#include "xform.h"
#include <opencv\highgui.h>

mat_s* create_init_img(mat_s*, int, double);
mat_s* convert_to_gray32(mat_s*);
mat_s*** build_gauss_pyr(mat_s*, int, int, double);
mat_s* downsample(mat_s*);
mat_s*** build_dog_pyr(mat_s***, int, int);
f_seq* scale_space_extrema(mat_s*** dog_pyr, int octvs, int intvls, double contr_thr, int curv_thr);
void calc_feature_scales(f_seq *features, double sigma, int intvls);
void adjust_for_img_dbl(f_seq *features);
mat_s* deriv_3D(mat_s*** dog_pyr, int octv, int intvl, int r, int c);
mat_s* hessian_3D(mat_s*** dog_pyr, int octv, int intvl, int r, int c);
void normalize_descr(struct feature* feat);
void interp_hist_entry(double*** hist, double rbin, double cbin, double obin, double mag, int d, int n);
void release_descr_hist(double**** hist, int d);
double*** descr_hist(mat_s *img, int r, int c, double ori, double scl, int d, int n);
double *ori_hist(mat_s *img, int r, int c, int n, int rad, double sigma);
int calc_grad_mag_ori(mat_s *img, int r, int c, double *mag, double *ori);
double dominant_ori(double *hist, int n);
void smooth_ori_hist(double *hist, int n);
void add_good_ori_features(f_seq *features, double *hist, int n, double mag_thr, struct feature* feat);
struct feature* clone_feature(struct feature* feat);
void hist_to_descr(double*** hist, int d, int n, struct feature* feat);
void calc_feature_oris(f_seq *features, mat_s*** gauss_pyr);
void compute_descriptors(f_seq* features, mat_s*** gauss_pyr, int d, int n);
int calc_sift_features(mat_s* img, struct feature** feat, int intvls, double sigma, double contr_thr, int curv_thr, int img_dbl, int descr_width, int descr_hist_bin);
void release_pyr(mat_s**** pyr, int octvs, int n);

void image_show_f32(char *str1, mat_s *mat, int one_sign)
{
	mat_s *mat_scale = NULL, *mat_u8 = NULL;
	if (one_sign)
	{
		mat_scale = mat_create(mat->rows, mat->cols, MAT_F32);
		mat_u8 = mat_create(mat->rows, mat->cols, MAT_U8);
		mat_scale = mat_convert_f32_scale(mat, 255.0f);
		mat_convert_f32_to_u8(mat_scale, mat_u8);
		image_show_mat(str1, mat_u8);
		mat_free(mat_scale);
		mat_free(mat_u8);
	}
	else
	{
		mat_u8 = mat_create(mat->rows, mat->cols, MAT_U8);
		mat_convert_f32_to_u8(mat, mat_u8);
		image_show_mat(str1, mat_u8);
		mat_free(mat_u8);
	}
}


int sift_features(mat_s *img,struct feature **feat)
{
	return calc_sift_features(img, feat, SIFT_INTVLS, SIFT_SIGMA, SIFT_CONTR_THR,
		SIFT_CURV_THR, SIFT_IMG_DBL, SIFT_DESCR_WIDTH,
		SIFT_DESCR_HIST_BINS);
}

int calc_sift_features(mat_s* img,struct feature** feat,int intvls,
						double sigma,double contr_thr,int curv_thr,
						int img_dbl,int descr_width,int descr_hist_bin)
{
	mat_s *init_img = NULL;
	mat_s ***gauss_pyr = NULL, ***dog_pyr = NULL;
	f_seq* features = NULL;
	int octvs;
	// step 1: build dog and dog_pyr
	init_img = create_init_img(img, img_dbl, sigma);
//	image_show_f32((char *)"init_img", init_img,1);
	octvs = (int)(log(MIN(init_img->rows, init_img->cols))/log(2) - 2);
	gauss_pyr = build_gauss_pyr(init_img, octvs, intvls, sigma);
	dog_pyr = build_dog_pyr(gauss_pyr, octvs, intvls);
	// step2: detect key_points
	features = scale_space_extrema(dog_pyr, octvs, intvls, contr_thr, curv_thr);
	calc_feature_scales(features,sigma,intvls);
	if (img_dbl)
	{
		adjust_for_img_dbl(features);
	}
	//step3: 
	calc_feature_oris(features, gauss_pyr);
	//step4: calc descriptors
	printf("step4 feature_seq length %d \n", feature_seq_length(features));
	compute_descriptors(features, gauss_pyr, descr_width, descr_hist_bin);
	int n = feature_seq_length(features);
	*feat = (feature*)calloc(n, sizeof(struct feature));
	*feat = feature_seq_to_array(features, *feat);
	for (int i = 0; i < n; i++)
	{
		free((*feat)[i].feature_data);
		(*feat)[i].feature_data = NULL;
	}
	feature_seq_clear(features);
	mat_free(init_img);
	release_pyr(&gauss_pyr, octvs, intvls + 3);
	release_pyr(&dog_pyr, octvs, intvls + 2);
	return n;
}


mat_s *create_init_img(mat_s *img,int img_dbl,double sigma)
{
	mat_s *gray, *dbl,*dbl_blur;
	float sig_diff;
	gray = convert_to_gray32(img);
	if (img_dbl)
	{
		sig_diff = (float)(sqrt(sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA * 4));
		dbl = mat_create(img->rows*2,img->cols*2,MAT_F32);
		mat_resize_double_f32(gray, dbl);
		dbl_blur = image_smooth_gaussian_rt(dbl,sig_diff,sig_diff);
		mat_free(gray);
		mat_free(dbl);
		return dbl_blur;
	}
	else
	{
		sig_diff = (float)(sqrt(sigma * sigma - SIFT_INIT_SIGMA * SIFT_INIT_SIGMA * 4));
		image_smooth_gaussian(gray, gray, sig_diff, sig_diff);
		return gray;
	}
}

mat_s *convert_to_gray32(mat_s *img)
{
	mat_s *gray8, *gray32;
	gray8 = mat_clone_u8_f32(img);
	gray32 = mat_convert_f32_scale(gray8,1.0f/255.0f);
	mat_free(gray8);
	return gray32;
}

mat_s ***build_gauss_pyr(mat_s *base, int octvs, int intvls, double sigma)
{
	mat_s ***gauss_pyr;
	double *sig = (double *)calloc(intvls + 3, sizeof(double));
	double sig_total, sig_prev, k;
	int i, o;
	gauss_pyr = (mat_s ***)calloc(octvs, sizeof(mat_s **));
	for (i = 0; i < octvs; i++)
	{
		gauss_pyr[i] = (mat_s**)calloc(intvls + 3, sizeof(mat_s*));
	}
	sig[0] = sigma;
	k = pow(2.0f, 1.0f / intvls);
	for (i = 1; i < intvls + 3; i++)
	{
		sig_prev = pow(k, i - 1) * sigma;
		sig_total = sig_prev * k;
		sig[i] = sqrt(sig_total * sig_total - sig_prev * sig_prev);
	}
	for (o = 0; o < octvs; o++)
	{
		for (i = 0; i < intvls + 3; i++)
		{
			if (o == 0 && i == 0)
			{
				gauss_pyr[o][i] = mat_clone_f32(base);
//				image_show_f32((char *)"o=0,i=0", gauss_pyr[o][i], 1);
			}
			else if (i == 0)
			{
				gauss_pyr[o][i] = mat_resize_half_f32(gauss_pyr[o-1][intvls]);
//				image_show_f32((char *)"i=0", gauss_pyr[o][i], 1);
			}
			else
			{
//				printf("o ====== %d i =======%d\n", o,i);
				gauss_pyr[o][i] = image_smooth_gaussian_rt(gauss_pyr[o][i-1],sig[i],sig[i]);
//				image_show_f32((char *)"o_i", gauss_pyr[o][i], 1);
			}
		}
	}
	free(sig);
	return gauss_pyr;
}

mat_s ***build_dog_pyr(mat_s ***gauss_pyr,int octvs,int intvls)
{
	mat_s ***dog_pyr;
	int i, o;
	dog_pyr = (mat_s ***)calloc(octvs, sizeof(mat_s**));
	for (i = 0; i < octvs; i++)
	{
		dog_pyr[i] = (mat_s**)calloc(intvls + 2, sizeof(mat_s*));
	}
	for (o = 0; o < octvs; o++)
	{
		for (i = 0; i < intvls + 2; i++)
		{
			dog_pyr[o][i] = mat_create(gauss_pyr[o][i]->rows, gauss_pyr[o][i]->cols, MAT_F32);
			mat_sub_r(gauss_pyr[o][i + 1], gauss_pyr[o][i], dog_pyr[o][i]);
//			printf("o ====== %d i =======%d\n", o, i);
//			image_show_f32((char *)"dog_o_i", dog_pyr[o][i], 1);
		}
	}
	return dog_pyr;
}

float mat_val(mat_s* mat, int row, int col)
{	
	float *mat_buf = (float *)mat->buf;
	return mat_buf[row*mat->cols + col];
}


int is_extremum(mat_s*** dog_pyr,int octv,int intvl,int r,int c)
{
	float val = mat_val(dog_pyr[octv][intvl], r, c);
	int i, j, k;
	if (val > 0)
	{
		for (i = -1; i <= 1; i++)
		{
			for (j = -1; j <= 1; j++)
			{
				for (k = -1; k <= 1; k++)
				{
					if (val < mat_val(dog_pyr[octv][intvl + i], r+j, c+k))
					{
						return 0;
					}
				}
			}
		}
	}
	else
	{
		for (i = -1; i <= 1; i++)
		{
			for (j = -1; j <= 1; j++)
			{
				for (k = -1; k <= 1; k++)
				{
					if (val > mat_val(dog_pyr[octv][intvl + i], r + j, c + k))
					{
						return 0;
					}
				}
			}
		}
	}
	return 1;
}


void interp_step(mat_s*** dog_pyr,int octv,int intvl,int r,int c,double* xi,double* xr,double* xc)
{
	mat_s* dD, *H, *H_inv, *X = NULL;
	dD = deriv_3D(dog_pyr, octv, intvl, r, c);
	H = hessian_3D(dog_pyr, octv, intvl, r, c);
	/*calc inv of H   use SVD method  */
	H_inv = svd_inv_f32(H);
	float *inv = mat_f32_to_array_f32(H_inv);
	float *d_D = mat_f32_to_array_f32(dD);
//	printf("inv \n");
//	mat_printf_f32(inv, H_inv->rows, H_inv->cols);
//	printf("d_D \n");
//	mat_printf_f32(d_D, dD->rows, dD->cols);
	*xc = -(inv[0] * d_D[0] + inv[1] * d_D[1] + inv[2] * d_D[2]);
	*xr = -(inv[3] * d_D[0] + inv[4] * d_D[1] + inv[5] * d_D[2]);
	*xi = -(inv[6] * d_D[0] + inv[7] * d_D[1] + inv[8] * d_D[2]);

//	X = mat_mul(H_inv,dD);
//	float* x_buf = (float *)X->buf;
//	*xi = -x_buf[2];
//	*xr = -x_buf[1];
//	*xc = -x_buf[0];
//	mat_printf(X);
//	printf("%f %f %f \n", -x_buf[2], -x_buf[1], -x_buf[0]);
//	printf("*xi = %f xr = %f xc = %f \n", *xi, *xr, *xc);

	mat_free(dD);
	mat_free(H);
	mat_free(H_inv);
	free(inv);
	free(d_D);
//	mat_free(X);


}

mat_s* deriv_3D(mat_s*** dog_pyr,int octv,int intvl,int r,int c)
{
	mat_s* dI;
	float dx, dy, ds;
	dx = (mat_val(dog_pyr[octv][intvl], r, c + 1) - mat_val(dog_pyr[octv][intvl], r, c - 1)) / 2.0f;
	dy = (mat_val(dog_pyr[octv][intvl], r + 1, c) - mat_val(dog_pyr[octv][intvl], r - 1, c)) / 2.0f;
	ds = (mat_val(dog_pyr[octv][intvl + 1], r, c) - mat_val(dog_pyr[octv][intvl - 1], r, c)) / 2.0f;
	dI = mat_create(3, 1, MAT_F32);
	mat_set_f32(dI, 0, 0, dx);
	mat_set_f32(dI, 1, 0, dy);
	mat_set_f32(dI, 2, 0, ds);
	return dI;
}

mat_s* hessian_3D(mat_s*** dog_pyr,int octv,int intvl,int r,int c)
{
	mat_s *H;
	float v, dxx, dyy, dss, dxy, dxs, dys;
	v = mat_val(dog_pyr[octv][intvl], r, c);
	dxx = (mat_val(dog_pyr[octv][intvl], r, c + 1) +
		mat_val(dog_pyr[octv][intvl], r, c - 1) - 2 * v);
	dyy = (mat_val(dog_pyr[octv][intvl], r + 1, c) +
		mat_val(dog_pyr[octv][intvl], r - 1, c) - 2 * v);
	dss = (mat_val(dog_pyr[octv][intvl + 1], r, c) +
		mat_val(dog_pyr[octv][intvl - 1], r, c) - 2 * v);
	dxy = (mat_val(dog_pyr[octv][intvl], r + 1, c + 1) -
		mat_val(dog_pyr[octv][intvl], r + 1, c - 1) +
		mat_val(dog_pyr[octv][intvl], r - 1, c + 1) +
		mat_val(dog_pyr[octv][intvl], r - 1, c - 1)) / 4.0f;
	dxs = (mat_val(dog_pyr[octv][intvl + 1], r, c + 1) -
		mat_val(dog_pyr[octv][intvl + 1], r, c - 1) -
		mat_val(dog_pyr[octv][intvl - 1], r, c + 1) +
		mat_val(dog_pyr[octv][intvl - 1], r, c - 1)) / 4.0f;
	dys = (mat_val(dog_pyr[octv][intvl + 1], r + 1, c) -
		mat_val(dog_pyr[octv][intvl + 1], r - 1, c) -
		mat_val(dog_pyr[octv][intvl - 1], r + 1, c) +
		mat_val(dog_pyr[octv][intvl - 1], r - 1, c)) / 4.0f;

	H = mat_create(3, 3, MAT_F32);
	mat_set_f32(H, 0, 0, dxx);
	mat_set_f32(H, 0, 1, dxy);
	mat_set_f32(H, 0, 2, dxs);
	mat_set_f32(H, 1, 0, dxy);
	mat_set_f32(H, 1, 1, dyy);
	mat_set_f32(H, 1, 2, dys);
	mat_set_f32(H, 2, 0, dxs);
	mat_set_f32(H, 2, 1, dys);
	mat_set_f32(H, 2, 2, dss);
	return H;
}

double interp_contr(mat_s*** dog_pyr, int octv, int intvl, int r,int c, double xi, double xr, double xc)
{
	mat_s* dD;
	double t = 0, x[3] = { xc,xr,xi };
	dD = deriv_3D(dog_pyr, octv, intvl, r, c);
	float *dD_buf = (float *)dD->buf;
	/* dD'*x =T (3*1)'*(3*1) = (1*1)  */
	t = dD_buf[0] * x[0] + dD_buf[1] * x[1] + dD_buf[2] * x[2];
	mat_free(dD);
	return mat_val(dog_pyr[octv][intvl], r, c) + 0.5f * t;
}

struct feature* new_feature()
{
	struct feature* feat;
	struct detection_data* ddata;
	feat = (feature*)malloc(sizeof(struct feature));
	memset(feat, 0, sizeof(struct feature));
	ddata = (detection_data*)malloc(sizeof(struct detection_data));
	memset(ddata, 0, sizeof(struct detection_data));
	feat->feature_data = ddata;
	feat->type = FEATURE_LOWE;
	return feat;
}


struct feature* interp_extremum(mat_s*** dog_pyr,int octv,int intvl,int r,int c,int intvls,double contr_thr)
{
	struct feature* feat;
	struct detection_data* ddata;
	double xi = 0, xr = 0, xc = 0, contr;
	int i = 0;
	while (i < SIFT_MAX_INTERP_STEPS)
	{
		interp_step(dog_pyr, octv, intvl, r, c, &xi, &xr, &xc);
		if (ABS(xi) < 0.5f && ABS(xr) < 0.5f && ABS(xc) < 0.5f)
		{
			break;
		}
		c += (int)(xc + 0.5f);
		r += (int)(xr + 0.5f);
		intvl += (int)(xi + 0.5f);
		if (intvl < 1 || intvl > intvls || c < SIFT_IMG_BORDER || r < SIFT_IMG_BORDER ||
			c >= dog_pyr[octv][0]->cols - SIFT_IMG_BORDER ||
			r >= dog_pyr[octv][0]->rows - SIFT_IMG_BORDER)
		{
			return NULL;
		}
		i++;
	}
	if (i >= SIFT_MAX_INTERP_STEPS)
	{
		return NULL;
	}
	contr = interp_contr(dog_pyr, octv, intvl, r, c, xi, xr, xc);
	if (ABS(contr) < contr_thr / intvls)
	{
		return NULL;
	}
	feat = new_feature();
	ddata = feat_detection_data(feat);
	feat->img_pt.x = feat->x = (c + xc) * pow(2.0f, octv);
	feat->img_pt.y = feat->y = (r + xr) * pow(2.0f, octv);
	ddata->r = r;
	ddata->c = c;
	ddata->octv = octv;
	ddata->intvl = intvl;
	ddata->sub_intvl = xi;
	return feat;

}

int is_too_edge_like(mat_s* dog_img,int r,int c,int curv_thr)
{
	double d, dxx, dyy, dxy, tr, det,tmp1,tmp2;
	d = mat_val(dog_img, r, c);
	dxx = mat_val(dog_img, r, c + 1) + mat_val(dog_img, r, c - 1) - 2 * d;
	dyy = mat_val(dog_img, r + 1, c) + mat_val(dog_img, r - 1, c) - 2 * d;
	dxy = (mat_val(dog_img, r + 1, c + 1) - mat_val(dog_img, r + 1, c - 1) -
		mat_val(dog_img, r - 1, c + 1) + mat_val(dog_img, r - 1, c - 1)) / 4.0f;
	tr = dxx + dyy;
	det = dxx * dyy - dxy * dxy;
	tmp1 = tr * tr / det;
	tmp2 = (curv_thr + 1.0f)*(curv_thr + 1.0f) / curv_thr;
	if (det <= 0)
	{
		return 1;
	}
	if (tmp1 < tmp2)
	{
		return 0;
	}
	return 1;

}

f_seq* scale_space_extrema(mat_s*** dog_pyr,int octvs,int intvls,double contr_thr,int curv_thr)
{
	f_seq* features = NULL;
	double prelim_contr_thr = 0.5f * contr_thr / intvls;
	struct feature* feat;
	struct detection_data* ddata;
	int o, i, r, c;
	for (o = 0; o < octvs; o++)
	{
		for (i = 1; i <= intvls; i++)
		{
			for (r = SIFT_IMG_BORDER; r < dog_pyr[o][0]->rows - SIFT_IMG_BORDER; r++)
			{
				for (c = SIFT_IMG_BORDER; c < dog_pyr[o][0]->cols - SIFT_IMG_BORDER; c++)
				{
					if (ABS(mat_val(dog_pyr[o][i], r, c)) > prelim_contr_thr)
					{
						if (is_extremum(dog_pyr, o, i, r, c))
						{
							feat = interp_extremum(dog_pyr, o, i, r, c, intvls, contr_thr);
							if (feat)
							{
	//							printf("o == %d i == %d r == %d c == %d \n", o, i, r, c);
								ddata = feat_detection_data(feat);
								if (!is_too_edge_like(dog_pyr[ddata->octv][ddata->intvl], ddata->r, ddata->c, curv_thr))
								{
									feature_seq_add_feature_1(&features,feat);
								}
								else
								{
									free(ddata);
								}
//								free(feat);
							}
						}
					}
				}
			}
		}
	}
	return features;
}

void calc_feature_scales(f_seq* features,double sigma,int intvls)
{
	struct feature* feat = NULL;
	struct detection_data* ddata;
	double intvl;
	int i, n = 0;
	n = feature_seq_length(features);
	for (i = 0;i < n; i++)
	{
		feat = feature_seq_get_specific_feature(features,i);
		ddata = feat_detection_data(feat);
		intvl = ddata->intvl + ddata->sub_intvl;
		feat->scl = sigma * pow(2.0f, ddata->octv + intvl / intvls);
		ddata->scl_octv = sigma * pow(2.0f, intvl / intvls);
	}
}

void adjust_for_img_dbl(f_seq *features)
{
	struct feature* feat = NULL;
	int i, n = 0;
	n = feature_seq_length(features);
	for (i = 0; i < n; i++)
	{
		feat = feature_seq_get_specific_feature(features,i);
		feat->x /= 2.0f;
		feat->y /= 2.0f;
		feat->scl /= 2.0f;
		feat->img_pt.x /= 2.0f;
		feat->img_pt.y /= 2.0f;
	}
}

void calc_feature_oris(f_seq *features,mat_s*** gauss_pyr)
{
	struct feature* feat;
	struct detection_data* ddata;
	double *hist;
	double omax;
	int i, j, n = 0;
	n = feature_seq_length(features);
	printf("feature_seq_length: %d \n", n);
	for (i = 0; i < n; i++)
	{
		feat = (struct feature*)malloc(sizeof(struct feature));
		feature_seq_pop_front(&features, feat);
//		printf("pop_front length %d\n", feature_seq_length(features));
		ddata = feat_detection_data(feat);
		hist = ori_hist(gauss_pyr[ddata->octv][ddata->intvl], ddata->r, ddata->c,SIFT_ORI_HIST_BINS,
			int(SIFT_ORI_RADIUS * ddata->scl_octv + 0.5f),SIFT_ORI_SIG_FCTR * ddata->scl_octv);
		for (j = 0; j < SIFT_ORI_SMOOTH_PASSES; j++)
		{
			smooth_ori_hist(hist, SIFT_ORI_HIST_BINS);
		}
		omax = dominant_ori(hist, SIFT_ORI_HIST_BINS);
		add_good_ori_features(features, hist, SIFT_ORI_HIST_BINS,omax * SIFT_ORI_PEAK_RATIO, feat);
//		printf("pop_front length %d\n", feature_seq_length(features));

		free(ddata);
		free(feat);
		free(hist);
	}
	n = feature_seq_length(features);
	printf("add good ori feature_seq_length: %d \n", n);


}

double *ori_hist(mat_s *img,int r,int c,int n,int rad,double sigma)
{
	double *hist;
	double mag, ori, w, exp_denom, PI2 = PI * 2.0f;
	int bin, i, j;
	hist = (double *)calloc(n, sizeof(double));
	exp_denom = 2.0f * sigma * sigma;
	for (i = -rad;i <= rad;i++)
	{
		for (j = -rad;j <= rad; j++)
		{
			if (calc_grad_mag_ori(img, r + i, c + j, &mag, &ori))
			{
				w = exp(-(i*i + j * j) / exp_denom);
				bin = (int)(n * (ori + PI) / PI2 + 0.5f);
				bin = (bin < n) ? bin : 0;
				hist[bin] += w * mag;
			}
		}
	}
	return hist;
}

int calc_grad_mag_ori(mat_s *img,int r,int c,double *mag,double *ori)
{
	double dx, dy;
	if (r > 0 && r < img->rows - 1 && c > 0 && c < img->cols - 1)
	{
		dx = mat_val(img, r, c + 1) - mat_val(img, r, c - 1);
		dy = mat_val(img, r - 1, c) - mat_val(img, r + 1, c);
		*mag = sqrt(dx * dx + dy * dy);
		*ori = atan2(dy, dx);
		return 1;
	}
	else
	{
		return 0;
	}
}

void smooth_ori_hist(double *hist,int n)
{
	double prev, tmp, h0 = hist[0];
	int i;
	prev = hist[n - 1];
	for (i = 0; i < n; i++)
	{
		tmp = hist[i];
		hist[i] = 0.25f * prev + 0.5f * hist[i] + 0.25f * ((i + 1 == n) ? h0 : hist[i + 1]);
		prev = tmp;
	}
}

double dominant_ori(double *hist, int n)
{
	double omax;
	int maxbin, i;
	omax = hist[0];
	maxbin = 0;
	for (i = 1; i < n; i++)
	{
		if (hist[i] > omax)
		{
			omax = hist[i];
			maxbin = i;
		}
	}
	return omax;
}

#define interp_hist_peak(l,c,r) (0.5f * (l-r)/(l-2.0f*c+r))

void add_good_ori_features(f_seq *features,double *hist,int n,double mag_thr,struct feature* feat)
{
	struct feature* new_feat;
	double bin, PI2 = PI * 2.0f;
	int l, r, i;
	for (i = 0; i < n; i++)
	{
		l = (i == 0) ? n - 1 : i - 1;
		r = (i + 1) % n;
		if (hist[i] > hist[l] && hist[i] > hist[r] && hist[i] >= mag_thr)
		{
			bin = i + interp_hist_peak(hist[l], hist[i], hist[r]);
			bin = (bin < 0) ? n + bin : (bin >= n) ? bin - n : bin;
			new_feat = clone_feature(feat);
			new_feat->ori = ((PI2 * bin) / n) - PI;
			feature_seq_add_feature_1(&features,new_feat);
//			free(new_feat);
		}
	}
}

struct feature* clone_feature(struct feature* feat)
{
	struct feature* new_feat;
	struct detection_data* ddata;
	new_feat = new_feature();
	ddata = feat_detection_data(new_feat);
	memcpy(new_feat, feat, sizeof(struct feature));
	memcpy(ddata, feat_detection_data(feat), sizeof(struct detection_data));
	new_feat->feature_data = ddata;
	return new_feat;
}

void compute_descriptors(f_seq* features,mat_s*** gauss_pyr,int d,int n)
{
	struct feature* feat = NULL;
	struct detection_data* ddata;
	double*** hist;
	int i, k = 0;
	k = feature_seq_length(features);
	for (i = 0; i < k; i++)
	{
		feat = feature_seq_get_specific_feature(features, i);
		ddata = feat_detection_data(feat);
		hist = descr_hist(gauss_pyr[ddata->octv][ddata->intvl], ddata->r, ddata->c,feat->ori, ddata->scl_octv, d, n);
		hist_to_descr(hist, d, n, feat);
		release_descr_hist(&hist, d);
	}
}

double*** descr_hist(mat_s *img, int r, int c, double ori,double scl, int d, int n)
{
	double*** hist;
	double cos_t, sin_t, hist_width, exp_denom, r_rot, c_rot, grad_mag;
	double grad_ori, w, rbin, cbin, obin, bins_per_rad, PI2 = 2.0f * PI;
	int radius, i, j;
	hist = (double***)calloc(d, sizeof(double **));
	for (i = 0; i < d; i++)
	{
		hist[i] = (double **)calloc(d, sizeof(double *));
		for (j = 0; j < d; j++)
		{
			hist[i][j] = (double *)calloc(n, sizeof(double));
		}
	}

	cos_t = cos(ori);
	sin_t = sin(ori);
	bins_per_rad = n / PI2;
	exp_denom = d * d * 0.5f;
	hist_width = SIFT_DESCR_SCL_FCTR * scl;
	radius = (int)(hist_width * sqrt(2) * (d + 1.0f) * 0.5f + 0.5f);

	for (i = -radius; i <= radius; i++)
	{
		for (j = -radius; j <= radius; j++)
		{
			c_rot = (j * cos_t - i * sin_t) / hist_width;
			r_rot = (j * sin_t + i * cos_t) / hist_width;
			rbin = r_rot + d / 2 - 0.5f;
			cbin = c_rot + d / 2 - 0.5f;
			if (rbin > -1.0f && rbin < d && cbin > -1.0f && cbin < d)
			{
				if (calc_grad_mag_ori(img, r + i, c + j, &grad_mag, &grad_ori))
				{
					grad_ori -= ori;
					while (grad_ori < 0.0f)
					{
						grad_ori += PI2;
					}
					while (grad_ori >= PI2)
					{
						grad_ori -= PI2;
					}
					obin = grad_ori * bins_per_rad;
					w = exp(-(c_rot * c_rot + r_rot * r_rot) / exp_denom);
					interp_hist_entry(hist, rbin, cbin, obin, grad_mag * w, d, n);
				}
			}
		}
	}
	return hist;
}

void interp_hist_entry(double*** hist,double rbin,double cbin,double obin,double mag,int d,int n)
{
	double d_r, d_c, d_o, v_r, v_c, v_o;
	double** row, *h;
	int r0, c0, o0, rb, cb, ob, r, c, o;
	r0 = (int)floor(rbin);
	c0 = (int)floor(cbin);
	o0 = (int)floor(obin);
	d_r = rbin - r0;
	d_c = cbin - c0;
	d_o = obin - o0;
	for (r = 0; r <= 1; r++)
	{
		rb = r0 + r;
		if (rb >= 0 && rb < d)
		{
			v_r = mag * ((r == 0) ? 1.0f - d_r : d_r);
			row = hist[rb];
			for (c = 0; c <= 1; c++)
			{
				cb = c0 + c;
				if (cb >= 0 && cb < d)
				{
					v_c = v_r * ((c == 0) ? 1.0f - d_c : d_c);
					h = row[cb];
					for (o = 0; o <= 1; o++)
					{
						ob = (o0 + o) % n;
						v_o = v_c * ((o == 0) ? 1.0f - d_o : d_o);
						h[ob] += v_o;
					}
				}
			}
		}
	}
}

void hist_to_descr(double*** hist, int d, int n, struct feature* feat)
{
	int int_val, i, r, c, o, k = 0;
	for (r = 0; r < d; r++)
	{
		for (c = 0; c < d; c++)
		{
			for (o = 0; o < n; o++)
			{
				feat->descr[k++] = hist[r][c][o];
			}
		}
	}
	feat->d = k;
	normalize_descr(feat);
	for (i = 0; i < k; i++)
	{
		if (feat->descr[i] > SIFT_DESCR_MAG_THR)
		{
			feat->descr[i] = SIFT_DESCR_MAG_THR;
		}
	}
	normalize_descr(feat);
	for (i = 0; i < k; i++)
	{
		int_val = (int)(SIFT_INT_DESCR_FCTR * feat->descr[i]);
		feat->descr[i] = MIN(255, int_val);
	}
}

void normalize_descr(struct feature* feat)
{
	double cur, len_inv, len_sq = 0.0f;
	int i, d = feat->d;
	for (i = 0; i < d; i++)
	{
		cur = feat->descr[i];
		len_sq += cur * cur;
	}
	len_inv = 1.0f / sqrt(len_sq);
	for (i = 0; i < d; i++)
	{
		feat->descr[i] *= len_inv;
	}
}

void release_descr_hist(double**** hist, int d)
{
	int i, j;
	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			free((*hist)[i][j]);
		}
		free((*hist)[i]);
	}
	free(*hist);
	*hist = NULL;
}

void release_pyr(mat_s**** pyr,int octvs,int n)
{
	int i, j;
	for (i = 0; i < octvs; i++)
	{
		for (j = 0; j < n;j++)
		{
			mat_free((*pyr)[i][j]);
		}
		free((*pyr)[i]);
	}
	free(*pyr);
	*pyr = NULL;
}


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
	mat = mat_from_data_u8(p, m, n);
//	image_show_mat(str,mat);
	free(p);
	return mat;
}


/* the maximum number of keypoint NN candidates to check during BBF search */
#define KDTREE_BBF_MAX_NN_CHKS 200

/* threshold on squared ratio of distances between NN and 2nd NN */
#define NN_SQ_DIST_RATIO_THR 0.49

IplImage* stacked_images(IplImage *img1, IplImage *img2)
{
	IplImage *stacked = cvCreateImage(cvSize(MAX(img1->width, img2->width),
		img1->height + img2->height), IPL_DEPTH_8U, 1);
	cvZero(stacked);
	cvSetImageROI(stacked, cvRect(0, 0, img1->width, img1->height));
	cvAdd(img1, stacked, stacked, NULL);
	cvSetImageROI(stacked, cvRect(0, img1->height, img2->width, img2->height));
	cvAdd(img2, stacked, stacked, NULL);
	cvResetImageROI(stacked);
	return stacked;
}


void test_sift_module()
{
	int m = 480;
	int n = 640;
	int num_matched = 0;
	int n1, n2;
	struct feature* feat1,*feat2,*feat;
	struct feature** nbrs;
	struct kd_node* kd_root;
	double d0, d1;
	CvPoint pt1, pt2;
	uint8_t *data1, *data2;
	mat_s *mat1,*mat2;
	data1 = read_data_u8_from_file_path((char *)"mat1.dat",m,n);
	data2 = read_data_u8_from_file_path((char *)"mat2.dat", n,m);
	mat1 = u8_data_to_mat_u8(data1, m, n);
	mat2 = u8_data_to_mat_u8(data2, n, m);
	n1 = sift_features(mat1, &feat1);
	n2 = sift_features(mat2, &feat2);
	IplImage *img1, *img2,*stacked;
	img1 = data_to_image(data1,m,n);
	img2 = data_to_image(data2,n,m);
	stacked = stacked_images(img1, img2);
	kd_root = kdtree_build(feat2, n2);
	for (int i = 0; i < n1; i++)
	{
		feat = feat1 + i;
		int k = kdtree_bbf_knn(kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS);
		if (k == 2)
		{
			d0 = descr_dist_sq(feat, nbrs[0]);
			d1 = descr_dist_sq(feat, nbrs[1]);
			if (d0 < d1 * NN_SQ_DIST_RATIO_THR)
			{
				pt1 = cvPoint(cvRound(feat->x), cvRound(feat->y));
				pt2 = cvPoint(cvRound(nbrs[0]->x), cvRound(nbrs[0]->y));
				pt2.y += img1->height;
				cvLine(stacked, pt1, pt2, CV_RGB(255, 0, 255), 1, 8, 0);
				num_matched++;
				feat1[i].fwd_match = nbrs[0];
			}
		}
		free(nbrs);
	}
	printf("Found %d total matches\n", num_matched);
	cvNamedWindow("Matches", 1);
	cvShowImage("Matches", stacked);
	cvWaitKey(0);
#if 1
	mat_s* H;
	CvMat *H_t;
	H = ransac_xform(feat1, n1, FEATURE_FWD_MATCH, 4, 0.01, 3.0f, NULL, NULL);
	mat_printf_f32((float *)H->buf, 3, 3);
	if (H)
	{
		H_t = mat_to_CvMat(H);
		for (int i = 0; i < H_t->height; i++)
		{
			float* pdata = (float *)(H_t->data.ptr + i * H_t->step);
			for (int j = 0; j < H_t->width; j++)
			{
				printf("%f \n", *pdata);
				*pdata++;
			}

		}
		IplImage *xformed;
		xformed = cvCreateImage(cvGetSize(img2), IPL_DEPTH_8U, 1);
		cvWarpPerspective(img1, xformed, H_t,CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS,cvScalarAll(0));
		cvNamedWindow("Xformed", 1);
		cvShowImage("Xformed", xformed);
		cvWaitKey(0);
		cvReleaseImage(&xformed);
		cvReleaseMat(&H_t);
	}
#endif



	cvReleaseImage(&stacked);
	free(data1);
	free(data2);
	mat_free(mat1);
	mat_free(mat2);
	kdtree_release(kd_root);
	free(feat1);
	free(feat2);


}
