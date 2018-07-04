#ifndef _XFORM_H_
#define _XFORM_H_
#include "mat.h"
#include "image_features.h"
#include <cxcore.h>

struct feature;
struct ransac_data
{
	void* orig_feat_data;
	int sampled;
};

#define RANSAC_ERR_TOL 3
#define RANSAC_INLIER_FRAC_EST 0.25
#define RANSAC_PROB_BAD_SUPP   0.10
#define feat_ransac_data(feat) ((struct ransac_data*)(feat)->feature_data);

mat_s* ransac_xform(struct feature* features, int n, int mtype, int m, double p_badxform, double err_tol, feature*** inliers, int *n_in);

#endif
