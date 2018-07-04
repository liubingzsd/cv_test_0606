#ifndef _MIN_PQ_H_
#define _MIN_PQ_H_
#include <stdlib.h>
#define MINPQ_INIT_NALLOCD  512

struct pq_node
{
	void *data;
	int key;
};

struct min_pq
{
	struct pq_node* pq_array;
	int nallocd;
	int n;
};

struct min_pq* minpq_init();
int minpq_insert(struct min_pq* min_pq, void* data, int key);
void *minpq_get_min(struct min_pq* min_pq);
void *minpq_extract_min(struct min_pq* min_pq);
void minpq_release(struct min_pq** min_pq);
#endif
