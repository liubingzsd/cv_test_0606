#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "svd_solver.h"
#include "xform.h"

static int get_matched_features(struct feature* features, int n, int mtype, struct feature*** matched);
static int calc_min_inliers(int n, int m, double p_badsupp, double p_badxform);
static struct feature** draw_ransac_sample(struct feature** features, int n, int m);
static void extract_corresp_pts(struct feature** features, int n, int mtype, point_2d_64f_s** pts, point_2d_64f_s** mpts);
static int find_consensus(struct feature** features, int n, int mtype, mat_s* M,double err_tol, struct feature*** consensus);
static __inline double log_factorial(int n);
static __inline struct feature* get_match(struct feature* feat, int mtype);
static __inline void release_mem(point_2d_64f_s* pts1, point_2d_64f_s* pts2, struct feature** features);
point_2d_64f_s persp_xform_pt(point_2d_64f_s pt, mat_s* T);
mat_s* least_square_homog(point_2d_64f_s* pts, point_2d_64f_s* mpts, int n);

mat_s* ransac_xform(struct feature* features, int n, int mtype,int m, double p_badxform, double err_tol,feature*** inliers, int *n_in)

{
	struct feature** matched, **sample, **consensus, **consensus_max = NULL;
	struct ransac_data* rdata;
	point_2d_64f_s *pts, *mpts;
	mat_s* M = NULL;
	double p, in_frac = RANSAC_INLIER_FRAC_EST;
	int i, nm, in, in_min, in_max = 0, k = 0;
	nm = get_matched_features(features, n, mtype, &matched);
	if (nm < m)
	{
		printf("ransac_xform() not enough matches to compute xform \n");
		goto end;
	}
	srand((unsigned int)time(NULL));
	in_min = calc_min_inliers(nm, m, RANSAC_PROB_BAD_SUPP, p_badxform);
	p = pow(1.0f - pow(in_frac, m), k);
	i = 0;
	while (p > p_badxform)
	{
		sample = draw_ransac_sample(matched, nm, m);
		extract_corresp_pts(sample, m, mtype, &pts, &mpts);
		M = least_square_homog(pts, mpts, m);
		if (!M)
		{
			goto iteration_end;
		}
		in = find_consensus(matched, nm, mtype, M, err_tol, &consensus);
		if (in > in_max)
		{
			if (consensus_max)
			{
				free(consensus_max);
			}
			consensus_max = consensus;
			in_max = in;
			in_frac = (double)in_max / nm;
		}
		else
		{
			free(consensus);
		}
		mat_free(M);
	
iteration_end:
		release_mem(pts, mpts, sample);
		p = pow(1.0f - pow(in_frac, m), ++k);
	}
	/*calc final transform based on best consensus set*/
	if (in_max >= in_min)
	{
		extract_corresp_pts(consensus_max, in_max, mtype, &pts, &mpts);
		M = least_square_homog(pts, mpts, in_max);
		in = find_consensus(matched, nm, mtype, M, err_tol, &consensus);
		mat_free(M);
		release_mem(pts, mpts, consensus_max);
		extract_corresp_pts(consensus, in, mtype, &pts, &mpts);
		M = least_square_homog(pts, mpts, in);
		if (inliers)
		{
			*inliers = consensus;
			consensus = NULL; 
		}
		if (n_in)
		{
			*n_in = in;
		}
		release_mem(pts, mpts, consensus);
	}
	else if (consensus_max)
	{
		if (inliers)
		{
			*inliers = NULL;
		}
		if (n_in)
		{
			*n_in = 0;
		}
		free(consensus_max);
	}

end:
	for (i = 0; i < nm; i++)
	{
		rdata = feat_ransac_data(matched[i]);
		matched[i]->feature_data = rdata->orig_feat_data;
		free(rdata);
	}
	free(matched);
	return M;
}

mat_s* direct_linear_transform_homog(point_2d_64f_s* pts, point_2d_64f_s* mpts, int n)
{
	mat_s* H, *A;
	double *D, *V,*h;
	int i;
	if (n < 4)
	{
		return NULL;
	}
	A = mat_create(2 * n, 9, MAT_F64);
	mat_set_zero_f64(A);
	for (i = 0; i < n; i++)
	{
		mat_set_f64(A, 2 * i, 3, -pts[i].x);
		mat_set_f64(A, 2 * i, 4, -pts[i].y);
		mat_set_f64(A, 2 * i, 5, -1.0);
		mat_set_f64(A, 2 * i, 6, mpts[i].y * pts[i].x);
		mat_set_f64(A, 2 * i, 7, mpts[i].y * pts[i].y);
		mat_set_f64(A, 2 * i, 8, mpts[i].y);
		mat_set_f64(A, 2 * i + 1, 0, pts[i].x);
		mat_set_f64(A, 2 * i + 1, 1, pts[i].y);
		mat_set_f64(A, 2 * i + 1, 2, 1.0);
		mat_set_f64(A, 2 * i + 1, 6, -mpts[i].x * pts[i].x);
		mat_set_f64(A, 2 * i + 1, 7, -mpts[i].x * pts[i].y);
		mat_set_f64(A, 2 * i + 1, 8, -mpts[i].x);
	}
	D = (double *)calloc(2*n * 2*n,sizeof(double));
	V = (double *)calloc(9 * 9, sizeof(double));
	mat_svd_decmp(A, D, V);   	
	h = mat_get_col(V, 9, 9, 8);
	H = array_f64_to_mat_f32(h,3,3);
	free(D);
	free(V);
	free(h);
	mat_free(A);

	return H;
}

mat_s* least_square_homog(point_2d_64f_s* pts, point_2d_64f_s* mpts, int n)
{
	mat_s *H, *A;
	double *B,*x;
	int i;
	if (n < 4)
	{
		printf("lsq_homog() too few points in lsq_homog \n");
		return NULL;
	}
	/* AX = B*/
	A = mat_create(2 * n, 8, MAT_F64);
	B = (double *)calloc(2 * n,sizeof(double));
	x = (double *)calloc(9, sizeof(double));
	mat_set_zero_f64(A);
	for (i = 0; i < n; i++)
	{
		mat_set_f64(A, i, 0, pts[i].x);
		mat_set_f64(A, i + n, 3, pts[i].x);
		mat_set_f64(A, i, 1, pts[i].y);
		mat_set_f64(A, i + n, 4, pts[i].y);
		mat_set_f64(A, i, 2, 1.0);
		mat_set_f64(A, i + n, 5, 1.0);
		mat_set_f64(A, i, 6, -pts[i].x * mpts[i].x);
		mat_set_f64(A, i, 7, -pts[i].y * mpts[i].x);
		mat_set_f64(A, i + n, 6, -pts[i].x * mpts[i].y);
		mat_set_f64(A, i + n, 7, -pts[i].y * mpts[i].y);
	}
	for (i = 0; i < n; i++)
	{
		B[i] = mpts[i].x;
		B[i + n] = mpts[i].y;
	}
	/*svd solve X*/
	mat_svd_solver_f64(A, B, x);
	x[8] = 1.0f;
	/*X trans to H*/
	H = array_f64_to_mat_f32(x, 3, 3);
	mat_printf_f64(x, 1, 9);
	mat_free(A);
	free(B);
	free(x);

	return H;
}

double dist_sq_2d(point_2d_64f_s p1, point_2d_64f_s p2)
{
	double x_diff = p1.x - p2.x;
	double y_diff = p1.y - p2.y;
	return x_diff * x_diff + y_diff * y_diff;
}

double homog_xfer_err(point_2d_64f_s pt,point_2d_64f_s mpt,mat_s* H)
{
	point_2d_64f_s xpt = persp_xform_pt(pt, H);
	return sqrt(dist_sq_2d(xpt, mpt));
}
/*
matrix T 3*3
*/
point_2d_64f_s persp_xform_pt(point_2d_64f_s pt, mat_s* T)
{	
	float* t_buf = (float *)T->buf;
	double xy[3] = { pt.x,pt.y,1.0f }, uv[3] = { 0 };
	point_2d_64f_s rslt;
	uv[0] = t_buf[0] * xy[0] + t_buf[1] * xy[1] + t_buf[2] * xy[2];
	uv[1] = t_buf[3] * xy[0] + t_buf[4] * xy[1] + t_buf[5] * xy[2];
	uv[2] = t_buf[6] * xy[0] + t_buf[7] * xy[1] + t_buf[8] * xy[2];
	rslt = point_2d_create(uv[0] / uv[2], uv[1] / uv[2]);
	return rslt;
}

static int get_matched_features(struct feature* features, int n, int mtype, struct feature*** matched)
{
	struct feature** _matched;
	struct ransac_data* rdata;
	int i, m = 0;
	_matched = (feature**)calloc(n, sizeof(struct feature*));
	for (i = 0; i < n; i++)
	{
		if (get_match(features + i, mtype))
		{
			rdata = (ransac_data *)malloc(sizeof(struct ransac_data));
			memset(rdata, 0, sizeof(struct ransac_data));
			rdata->orig_feat_data = features[i].feature_data;
			_matched[m] = features + i;
			_matched[m]->feature_data = rdata;
			m++;
		}
		*matched = _matched;
	}
	return m;
}

static __inline struct feature* get_match(struct feature* feat, int mtype)
{
	if (mtype == FEATURE_MDL_MATCH)
	{
		return feat->mdl_match;
	}
	if (mtype == FEATURE_BCK_MATCH)
	{
		return feat->bck_match;
	}
	if (mtype == FEATURE_FWD_MATCH)
	{
		return feat->fwd_match;
	}
	return NULL;
}

static int calc_min_inliers(int n, int m, double p_badsupp, double p_badxform)
{
	double pi, sum;
	int i, j;
	for (j = m + 1; j <= n; j++)
	{
		sum = 0;
		for (i = j; i <= n; i++)
		{
			pi = (i - m)*log(p_badsupp) + (n - i + m)*log(1.0f - p_badsupp) +
				log_factorial(n - m) - log_factorial(i - m) - log_factorial(n - i);
			sum += exp(pi);
		}
		if (sum < p_badxform)
		{
			break;
		}
	}
	return j;
}

static __inline double log_factorial(int n)
{
	double f = 0;
	int i;
	for (i = 1; i <= n; i++)
	{
		f += log(i);
	}
	return f;
}

static struct feature** draw_ransac_sample(struct feature** features, int n, int m)
{
	struct feature** sample, *feat;
	struct ransac_data* rdata;
	int i, x;
	for (i = 0; i < n; i++)
	{
		rdata = feat_ransac_data(features[i]);
		rdata->sampled = 0;
	}
	sample = (feature **)calloc(m, sizeof(struct feature*));
	for (i = 0; i < m; i++)
	{
		do
		{
			x = rand() % n;
			feat = features[x];
			rdata = feat_ransac_data(feat);
		} while (rdata->sampled);
		sample[i] = feat;
		rdata->sampled = 1;
	}
	return sample;
}

static void extract_corresp_pts(struct feature** features, int n, int mtype,point_2d_64f_s** pts, point_2d_64f_s** mpts)
{
	struct feature* match;
	point_2d_64f_s* _pts, *_mpts;
	int i;
	_pts = (point_2d_64f_s *)calloc(n, sizeof(point_2d_64f_s));
	_mpts = (point_2d_64f_s *)calloc(n, sizeof(point_2d_64f_s));
	if (mtype == FEATURE_MDL_MATCH)
	{
		for (i = 0; i < n; i++)
		{
			match = get_match(features[i], mtype);
			if (!match)
			{
				printf("feature does not have match of type");
			}
			_pts[i] = features[i]->img_pt;
			_mpts[i] = match->mdl_pt;
		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			match = get_match(features[i], mtype);
			if (!match)
			{
				printf("feature does not have match of type");
			}
			_pts[i] = features[i]->img_pt;
			_mpts[i] = match->img_pt;
		}
	}
	*pts = _pts;
	*mpts = _mpts;
}

static int find_consensus(struct feature** features, int n, int mtype, mat_s* M,double err_tol, struct feature*** consensus)
{
	struct feature** _consensus;
	struct feature* match;
	point_2d_64f_s pt, mpt;
	double err;
	int i, in = 0;
	_consensus = (feature**)calloc(n, sizeof(struct feature*));
	if (mtype == FEATURE_MDL_MATCH)
	{
		for (i = 0; i < n; i++)
		{
			match = get_match(features[i], mtype);
			if (!match)
			{
				printf("feature does not have match of type \n");
			}
			pt = features[i]->img_pt;
			mpt = match->mdl_pt;
			err = homog_xfer_err(pt, mpt, M);
			if (err <= err_tol)
			{
				_consensus[in++] = features[i];
			}
		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			match = get_match(features[i], mtype);
			if (!match)
			{
				printf("feature does not have match of type \n");
			}
			pt = features[i]->img_pt;
			mpt = match->img_pt;
			err = homog_xfer_err(pt, mpt, M);
			if (err <= err_tol)
			{
				_consensus[in++] = features[i];
			}
		}
	}
	*consensus = _consensus;
	return in;
}

static __inline void release_mem(point_2d_64f_s* pts1, point_2d_64f_s* pts2, struct feature** features)
{
	free(pts1);
	free(pts2);
	if (features)
	{
		free(features);
	}
}