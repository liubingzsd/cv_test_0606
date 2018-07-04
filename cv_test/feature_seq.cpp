#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "feature_seq.h"

f_seq* feature_seq_add_feature(f_seq* feature_head,struct feature* feature)
{
	f_seq* temp, *list;
	list = (f_seq*)malloc(sizeof(f_seq));
	list->feat = feature;
	list->next = NULL;
	if (feature_head == NULL)
	{
		feature_head = list;
		feature_head->next = NULL;
	}
	else
	{
		temp = feature_head;
		while (temp)
		{
			if (temp->next == NULL)
			{
				temp->next = list;
				temp->next->next = NULL;
			}
			temp = temp->next;
		}
	}
	return feature_head;
}

void feature_seq_add_feature_1(f_seq** feature_head, struct feature* feature)
{
	f_seq* temp, *list;
	list = (f_seq*)malloc(sizeof(f_seq));
	list->feat = feature;
	list->next = NULL;
	if (*feature_head == NULL)
	{
		*feature_head = list;
		(*feature_head)->next = NULL;
	}
	else
	{
		temp = *feature_head;
		while (temp)
		{
			if (temp->next == NULL)
			{
				temp->next = list;
				temp->next->next = NULL;
			}
			temp = temp->next;
		}
	}
}

struct feature* feature_seq_pop_feature(f_seq** seq)
{
	struct feature* feat;
	f_seq* temp;
	feat = (*seq)->feat;
	temp = (*seq)->next;
	free(*seq);
	*seq = NULL;
	*seq = temp;
	return feat;
}

void feature_seq_pop_front(f_seq** seq,struct feature* feat)
{
	f_seq* temp = NULL;
	memcpy(feat, (*seq)->feat, sizeof(struct feature));
	temp = (*seq)->next;
	free(*seq);
	*seq = NULL;
	*seq = temp;
}

struct feature* feature_seq_get_specific_feature(f_seq* seq, int i)
{
	int num = 0;
	f_seq* temp = NULL;
	temp = seq;
	while (temp != NULL)
	{
		if (i == num)
		{
			break;
		}
		else
		{
			temp = temp->next;
		}
		num++;
	}
	return temp->feat;
}


int feature_seq_length(f_seq* seq)
{
	int num = 0;
	f_seq* temp;
	temp = seq;
	while (temp != NULL)
	{
		num++;
		temp = temp->next;
	}
//	printf("list length is %d \n", num);
	return num;
}

void feature_seq_clear(f_seq* seq)
{
	f_seq* temp;
	int n = 0;
	while (seq != NULL)
	{
		n++;
		temp = seq->next;
		free(seq);
		seq = NULL;
		seq = temp;
	}
	printf("clear %d element \n", n);
}

feature* feature_seq_to_array(f_seq* features,feature* feat)
{
	int i, n;
	feature* tmp;
	n = feature_seq_length(features);
	for (i = 0; i < n; i++)
	{
		tmp = feature_seq_get_specific_feature(features, i);
		memcpy(&feat[i], tmp, sizeof(struct feature));
	}
	return feat;
}


void test_feature_seq()
{
	int i;
	f_seq *feature_head = NULL;
	printf("__________test_feature_seq_add_feature___________\n");
	for (i = 0; i < 20; i++)
	{
		printf("i = %d \n", i);
		struct feature* feat = (struct feature *)malloc(sizeof(struct feature));
		feat->x = i;
		feat->y = i + 20;
//		printf("feat addr = 0x%x\n", feat);	
		feature_head = feature_seq_add_feature(feature_head,feat);
		printf("feat->x = %f feat->f = %f \n", feat->x, feat->y);
	}
	printf("__________test_feature_seq_get_specific_feature___________\n");
	for (i = 0; i < 20; i++)
	{
		printf("i = %d \n", i);
		struct feature* feat;
		feat = feature_seq_get_specific_feature(feature_head, i);
//		printf("feat addr = 0x%x\n", feat);
		printf("feat->x = %f feat->f = %d \n", feat->x, feat->y);
	}

	printf("__________test_feature_seq_to_array___________\n");
	feature* feat = (feature*)calloc(20, sizeof(struct feature));
	feat = feature_seq_to_array(feature_head,feat);
	for (i = 0; i < 20; i++)
	{
		printf("%f %f \n", feat[i].x, feat[i].y);
	}

	printf("__________test_feature_seq_pop_front___________\n");
	for (i = 0; i < 20; i++)
	{
		printf("i = %d \n", i);
		struct feature* feat;
		feat = (struct feature*)malloc(sizeof(struct feature));
		feature_seq_pop_front(&feature_head,feat);
		printf("feat->f = %f feat->f = %f \n", feat->x, feat->y);
		printf("feature length %d \n", feature_seq_length(feature_head));
	}

//	feature_seq_clear(feature_head);
//	printf("pop feature \n");
//	for (i = 0; i < 20; i++)
//	{
//		struct feature* feat ;
//		feat = feature_seq_pop_feature(&feature_head);
//		printf("feat addr = 0x%x\n", feat);
//	}

}