#ifndef _IMAGE_FEATURES_H_
#define _IMAGE_FEATURES_H_

/** FEATURE_OXFD <BR> FEATURE_LOWE */
enum feature_type
{
	FEATURE_OXFD,
	FEATURE_LOWE,
};

/** FEATURE_FWD_MATCH <BR> FEATURE_BCK_MATCH <BR> FEATURE_MDL_MATCH */
enum feature_match_type
{
	FEATURE_FWD_MATCH,
	FEATURE_BCK_MATCH,
	FEATURE_MDL_MATCH,
};


#define FEATURE_MAX_D 128

typedef struct 
{
	double x;
	double y;
}point_2d_64f_s;


struct feature
{
	double x;							/**< x coord */
	double y;							/**< y coord */
	double a;							/**< Oxford-type affine region parameter */
	double b;							/**< Oxford-type affine region parameter */
	double c;							/**< Oxford-type affine region parameter */
	double scl;							/**< scale of a Lowe-style feature */
	double ori;							/**< orientation of a Lowe-style feature */
	int d;								/**< descriptor length */
	double descr[FEATURE_MAX_D];		/**< descriptor */
	int type;							/**< feature type, OXFD or LOWE */
	int category;						/**< all-purpose feature category */
	struct feature* fwd_match;			/**< matching feature from forward image */
	struct feature* bck_match;			/**< matching feature from backmward image */
	struct feature* mdl_match;			/**< matching feature from model */
	point_2d_64f_s img_pt;				/**< location in image */
	point_2d_64f_s mdl_pt;				/**< location in model */
	void* feature_data;					/**< user-definable data */
};

double descr_dist_sq(struct feature* f1, struct feature* f2);
point_2d_64f_s point_2d_create(double x, double y);

#endif