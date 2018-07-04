#include <stdio.h>
#include <stdlib.h>
#include "linked_list.h"

struct list {
	int id;
	struct list* next;
};
static struct list* list_head = NULL;
static int id = 0;

static void list_add(struct list** head,struct list* list)
{
	struct list* temp;
	if (NULL == *head)
	{
		*head = list;
		(*head)->next = NULL;
	}
	else
	{
		temp = *head;
		while (temp)
		{
			if (NULL == temp->next)
			{
				temp->next = list;
				list->next = NULL;
			}
			temp = temp->next;
		}
	}
}

static void list_print(struct list **head)
{
	struct list *temp;
	temp = *head;
	printf("list information : \n");
	while (temp)
	{
//		printf("\t list %d : %x \n", temp->id, temp->next);
		temp = temp->next;
	}
}
static int list_delete(struct list** head,int id)
{
	struct list* temp, *p;
	temp = *head;
	if (NULL == temp)
	{
		printf("list is empty \n");
		return -1;
	}
	else
	{
		if (id == temp->id)
		{
			*head = temp->next;
			free(temp);
			return 0;
		}
		else
		{
			while (temp->next)
			{
				p = temp;
				temp = temp->next;
				if (id == temp->id)
				{
					p->next = temp->next;
					free(temp);
					return 0;
				}
			}
			return -1;
		}
	}
	return -1;
}

static int list_length(struct list** head)
{
	int num = 0;
	struct list* temp, *p;
	temp = *head;
	while (temp != NULL)
	{
		num++;
		p = temp->next;
		temp = p;
	}
	return num;
}

static void list_pop(struct list** head,struct list* pop)
{
	struct list* temp;
	pop->id = (*head)->id;
	pop->next = (*head)->next;
	temp = (*head)->next;
	free(*head);
	*head = NULL;
	*head = temp;
}
void test_pop(struct list** head)
{
	struct list* pop = NULL;
	pop = (struct list*)malloc(sizeof(struct list));
	list_pop(head, pop);
	printf("pop id = %d \n", pop->id);
}






int test_list()
{
	int i = 0;
	struct list* lists = NULL;
	struct list* pop = NULL;
	for (i = 0; i < 10; i++)
	{
		lists = (struct list*)malloc(sizeof(struct list));
		if (lists == NULL)
		{
			printf("malloc error \n");
		}
//		printf("list_id = %d \n", id);
		lists->id = id++;
//		printf("lists->id = %d \n", lists->id);
		list_add(&list_head, lists);
	}
	for (i = 0; i < 10; i++)
	{
		test_pop(&list_head);
	}
//	list_delete(&list_head, 3);
//	list_print(&list_head);
	return 0;
}
