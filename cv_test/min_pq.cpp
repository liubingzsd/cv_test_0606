#include <stdio.h>
#include "min_pq.h"

int array_double(void **array, int n, int size);
static void decrease_pq_node_key(struct pq_node* pq_array, int i, int key);
static void restore_minpq_order(struct pq_node* pq_array, int i, int n);

static __inline int parent(int i)
{
	return (i - 1) / 2;
}

static __inline int right(int i)
{
	return 2 * i + 2;
}

static __inline int left(int i)
{
	return 2 * i + 1;
}


struct min_pq* minpq_init()
{
	struct min_pq* min_pq;
	min_pq = (struct min_pq *)malloc(sizeof(struct min_pq));
	min_pq->pq_array = (pq_node *)calloc(MINPQ_INIT_NALLOCD,sizeof(struct pq_node));
	min_pq->nallocd = MINPQ_INIT_NALLOCD;
	min_pq->n = 0;
	return min_pq;
}

int minpq_insert(struct min_pq* min_pq, void *data, int key)
{
	int n = min_pq->n;

	if (min_pq->nallocd == n)
	{
		min_pq->nallocd = array_double((void **)&min_pq->pq_array, min_pq->nallocd, sizeof(struct pq_node));
		if (!min_pq->nallocd)
		{
			printf("__minpq_insert() unable to allocate memory \n");
			return 1;
		}
	}
	min_pq->pq_array[n].data = data;
	min_pq->pq_array[n].key = INT_MAX;
	decrease_pq_node_key(min_pq->pq_array, min_pq->n, key);
	min_pq->n++;
	return 0;
}

int array_double(void **array, int n, int size)
{
	void* tmp;
	tmp = realloc(*array, 2 * n * size);
	if (!tmp)
	{
		printf("array_double() unable to allocate memory \n");
		if (*array)
		{
			free(*array);
			*array = NULL;
		}
		return 0;
	}
	*array = tmp;
	return n * 2;
}

static void decrease_pq_node_key(struct pq_node* pq_array,int i,int key)
{
	struct pq_node tmp;
	if (key > pq_array[i].key)
	{
		return;
	}
	pq_array[i].key = key;
	while (i > 0 && pq_array[i].key < pq_array[parent(i)].key)
	{
		tmp = pq_array[parent(i)];
		pq_array[parent(i)] = pq_array[i];
		pq_array[i] = tmp;
		i = parent(i);
	}
}

void* minpq_extract_min(struct min_pq* min_pq)
{
	void *data;
	if (min_pq->n < 1)
	{
		printf("minpq_extract_min() pq empty \n");
		return NULL;
	}
	data = min_pq->pq_array[0].data;
	min_pq->n--;
	min_pq->pq_array[0] = min_pq->pq_array[min_pq->n];
	restore_minpq_order(min_pq->pq_array, 0, min_pq->n);
	return data;
}

static void restore_minpq_order(struct pq_node* pq_array, int i, int n)
{
	struct pq_node tmp;
	int l, r, min = i;
	l = left(i);
	r = right(i);
	if (l < n)
	{
		if (pq_array[l].key < pq_array[i].key)
		{
			min = l;
		}
	}
	if (r < n)
	{
		if (pq_array[r].key < pq_array[min].key)
		{
			min = r;
		}
	}
	if (min != i)
	{
		tmp = pq_array[min];
		pq_array[min] = pq_array[i];
		pq_array[i] = tmp;
		restore_minpq_order(pq_array, min, n);
	}
}

void minpq_release(struct min_pq** min_pq)
{
	if (!min_pq)
	{
		printf("minpq_release() NULL pointer error \n");
		return;
	}
	if (*min_pq && (*min_pq)->pq_array)
	{
		free((*min_pq)->pq_array);
		free(*min_pq);
		*min_pq = NULL;
	}
}