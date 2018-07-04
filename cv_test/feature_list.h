#ifndef _FEATURE_LIST_H_
#define _FEATURE_LIST_H_
#include "image_features.h"

struct feature_list
{
	struct feature* feat;
	struct feature_list* next;
};
typedef feature_list f_list;

 void feature_list_add_list(f_list** f_head, f_list* list);
 void feature_list_add_feature(f_list** f_head, struct feature* feat);
 int feature_list_length(f_list** f_head);
 void feature_list_clear(f_list** f_head);
 void feature_list_pop_feature(f_list** f_head, struct feature** feat);
 void test_f_list();
 f_list* feature_list_add_feat(struct feature* feat);

#endif