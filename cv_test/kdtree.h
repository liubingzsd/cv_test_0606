#ifndef _KDTREE_H_
#define _KDTREE_H_

struct feature;
struct kd_node
{
	int ki;
	double kv;
	int leaf;
	struct feature* features;
	int n;
	struct kd_node* kd_left;
	struct kd_node* kd_right;
};

struct kd_node* kdtree_build(struct feature* features, int n);
int kdtree_bbf_knn(struct kd_node* kd_root, struct feature* feat, int k,
	struct feature*** nbrs, int max_nn_chks);

void kdtree_release(struct kd_node* kd_root);
#endif