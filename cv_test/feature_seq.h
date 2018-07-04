#ifndef _FEATURE_SEQ_H_
#define _FEATURE_SEQ_H_
#include "image_features.h"

struct feature_seq
{
	struct feature* feat;
	struct feature_seq* next;
};
typedef feature_seq f_seq;

f_seq* feature_seq_add_feature(f_seq* seq,struct feature* feature);
void feature_seq_add_feature_1(f_seq** feature_head, struct feature* feature);
feature* feature_seq_to_array(f_seq* features, feature* feat);
void feature_seq_pop_front(f_seq** seq, struct feature* feat);
struct feature* feature_seq_pop_feature(f_seq** seq);
struct feature* feature_seq_get_specific_feature(f_seq* seq, int i);
int feature_seq_length(f_seq* seq);
void feature_seq_clear(f_seq* seq);
void test_feature_seq();

#endif
