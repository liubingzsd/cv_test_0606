#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "feature_list.h"

f_list* f_head = NULL;
f_list *feature_head = NULL;

void feature_list_add_list(f_list** f_head, f_list* list)
{
	f_list* temp;
	if (*f_head == NULL)
	{
		*f_head = list;
		(*f_head)->next = NULL;
	}
	else
	{
		temp = *f_head;
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

void feature_list_add_feature(f_list** f_head, struct feature* feat)
{
	static f_list* temp;
	if (*f_head == NULL)
	{
		*f_head = (f_list*)malloc(sizeof(f_list));
		(*f_head)->feat = feat;
		(*f_head)->next = NULL;
		temp = *f_head;
	}
	else
	{
		if (temp->next == NULL)
		{
			temp->next = (f_list*)malloc(sizeof(f_list));
			temp->next->feat = feat;
			temp->next->next = NULL;
		}
		temp = temp->next;
	}
}


void feature_list_add_feature_test(f_list** f_head, struct feature* feat)
{
	f_list* temp,*list;
	list = (f_list*)malloc(sizeof(f_list));
	list->feat = feat;
	list->next = NULL;
	if (*f_head == NULL)
	{
		*f_head = list;
		(*f_head)->next = NULL;
	}
	else
	{
		temp = *f_head;
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

f_list* feature_list_add_feat(struct feature* feat)
{
	f_list* temp, *list;
	list = (f_list*)malloc(sizeof(f_list));
	list->feat = feat;
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



int feature_list_length(f_list** f_head)
{
	int num = 0;
	f_list* temp, *p;
	temp = *f_head;
	while (temp != NULL)
	{
		num++;
		p = temp->next;
		temp = p;
	}
	return num;
}

int feature_list_length_1(f_list* f_head)
{
	int num = 0;
	f_list* temp, *p;
	temp = f_head;
	while (temp != NULL)
	{
		num++;
		p = temp->next;
		temp = p;
	}
	printf("list length is %d \n", num);
	return num;
}


void feature_list_pop_feature(f_list** f_head, struct feature** feat)
{
	f_list* temp;
	*feat = (*f_head)->feat;
	temp = (*f_head)->next;
	free(*f_head);
	*f_head = NULL;
	*f_head = temp;
}

struct feature* feature_list_pop_feature_1(f_list** feature_head)
{
	struct feature* feat;
	f_list* temp;
	feat = (*feature_head)->feat;
	temp = (*feature_head)->next;
	free(*feature_head);
	*feature_head = NULL;
	*feature_head = temp;
	return feat;
}



void feature_list_clear(f_list** f_head)
{
	f_list* temp;
	int n = 0;
	while(*f_head != NULL)
	{
		n++;
		temp = (*f_head)->next;
		free(*f_head);
		*f_head = NULL;
		*f_head = temp;
	}
	printf("clear %d element \n", n);
}

void feature_list_clear_1(f_list* f_head)
{
	f_list* temp;
	int n = 0;
	while (f_head != NULL)
	{
		n++;
		temp = f_head->next;
		free(f_head);
		f_head = NULL;
		f_head = temp;
	}
	printf("clear %d element \n", n);
}

void feature_list_print(f_list** f_head)
{
	f_list* temp;
	temp = *f_head;
	while (temp)
	{
//		printf("f_list feat = 0x%x f_list next = 0x%x \n", temp->feat, temp->next);
		temp = temp->next;
	}
}

void test_list_pop(f_list* feature_head)
{
	int i;
	int n = feature_list_length(&feature_head);
	for (i = 0; i < n; i++)
	{
		struct feature *feat1 = NULL;
		feature_list_pop_feature(&feature_head, &feat1);
//		printf("feat = 0x%x \n", feat1);
	}
}

void test_list_pop_1(f_list* feature_head)
{
	int i;
	int n = feature_list_length_1(feature_head);
	for (i = 0; i < n; i++)
	{
		struct feature *feat1 = NULL;
		feat1 = feature_list_pop_feature_1(&feature_head);
//		printf("feat = 0x%x \n", feat1);
	}
}


void test_f_list()
{
	int i;

	for (i = 0; i < 20; i++)
	{
		printf("i = %d \n", i);
		struct feature* feat = (struct feature *)malloc(sizeof(struct feature));
//		printf("feat addr = 0x%x\n", feat);
		feature_head = feature_list_add_feat(feat);
	}
	feature_list_print(&feature_head);
	feature_list_clear_1(feature_head);
//	test_list_pop_1(feature_head); 

}
