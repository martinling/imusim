/******************************************************************************
* utils.c
*
* Written by Ross Hemsley for McStas, September 2009.
* All code below is required for the intperolation functions, though many of 
* of the routines may be considered general enough for use anywhere they are
* required. Current utils include an array based list implementation, a doubly
* linked list implementation and an array based stack implementation.
*
*******************************************************************************/
#include "utils.h"

/* Unit Testing */
// #define _TEST_

/*******************************************************************************
* Array List implementation. We use this whenever we want arrays which can be
* resized dynamically. We expect O(1) amortised insert and extraction from this
* implementation. Other opertations are O(n).
*******************************************************************************/

int addToArrayList(arrayList *l, void* element)
{
   if (l->num_elements >= l->num_slots)
   {
   
      // we have to allocate more space
      l->num_slots *= 2;
      l->arr = realloc(l->arr, (l->num_slots*sizeof(void*)));
      // check that we haven't run out of memory
      if (l->arr == NULL)
      {
         fprintf(stderr, "Error: Out of Memory.\n");
         return -1;
      }
   }
   // add the element
   l->arr[l->num_elements] = element;
   l->num_elements++;
   
   // return the index where this element can be found.
   return (l->num_elements -1);
}

/******************************************************************************/

arrayList *newArrayList()
{
   arrayList *l;
   l = malloc(sizeof(arrayList));
   l->num_elements = 0;
   l->num_slots = 16;
   l->arr = malloc(16*sizeof(void*));
   return l;
}

/******************************************************************************/
// If this element is in the list, it will return the index of that element.
// Otherwise we return -1.

int arrayListGetIndex(arrayList *l, void *e)
{
  int i;
  for (i=0; i<arrayListSize(l); i++)
    if (getFromArrayList(l, i) == e) return i;
  return -1;
}

/******************************************************************************/

void** getArrayFromArrayList(arrayList *l)
{
   return l->arr;
}

/******************************************************************************/

// a special function, which only works when the arrayList contains only int*
int arrayListContains(arrayList * l , void * e)
{
   int i;
   for (i=0; i< l->num_elements; i++)
      if (e == l->arr[i]) return 1;
   return 0;
}

/******************************************************************************/

int arrayListSize(arrayList *l)
{
   return l->num_elements;
}

/******************************************************************************/

void * getFromArrayList (arrayList *l, int index)
{
   if(index >= 0 && index <  l->num_elements)
      return l->arr[index];
      
   return NULL;
}

/******************************************************************************/

void freeElements(arrayList *l)
{
  int i;
  
  for (i=0; i<arrayListSize(l); i++)
    free(getFromArrayList(l,i));

}

/******************************************************************************/
// We keep the memory associated with this list, but set it's number of elements
// to zero. This is effectively the same function as a memory heap.

void emptyArrayList(arrayList *l)
{
  l->num_elements=0; 
}

/******************************************************************************/

void freeArrayList(arrayList *l, void (*destructor)(void *e))
{
  int i;
  
  if (destructor)
    for (i=0;i<arrayListSize(l); i++)
      destructor(getFromArrayList(l,i));
  
  free(l->arr);
  free(l);
}

/******************************************************************************
* Doubly-linked list implementation. We use this implementation of a list when
* we want to have faster extractions from very large lists. Using this
* implementation we expect O(1) from all operations, except accessing an element
* by index, which is O(n). In most cases this will be much slower than using
* the array list implementation. (Except fast point removal in large sets) 
* because this implementation will not take advantage of the cache like the
* array list does.
*******************************************************************************/

linkedList *newLinkedList()
{
  linkedList *l = malloc(sizeof(linkedList));
    
  l->deadNodes = newStack();
  l->head  = NULL;
  l->last  = NULL;
  l->nelem = 0;
  
  return l;
}

/******************************************************************************/

listNode* addToLinkedList(linkedList *l, void *e)
{

  listNode *ln;
  
  if (!isEmpty(l->deadNodes))
    ln = pop(l->deadNodes);
  else
    ln = malloc(sizeof(listNode));

  ln->data = e;
  ln->next = NULL;
  ln->prev = l->last;
  
  // Put this element on the end. If this is the first element
  // then set the head.
  if (l->head)
    l->last->next = ln;
  else
    l->head = ln;
  l->last = ln;
  
  l->nelem ++;
  return ln;
}

/******************************************************************************/

void *getFromLinkedList(linkedList *l, int x)
{

  #ifdef DEBUG
  printf("Note: use of access by index in doubly linked list. Is this really "
         "what you want?\n");
  #endif
  
  if (! (0 <= x && x < linkedListSize(l)) )
  {
    fprintf(stderr, "list index out of bounds, linkedList-size: %d, index: %d.\n", 
      linkedListSize(l), x);
    exit(1);
  }
  
  listNode *thisNode = topOfLinkedList(l);
  int i;
  for (i=0; i<x; i++)
    thisNode = thisNode->next;
  return thisNode->data;

}

/******************************************************************************/

int linkedListSize(linkedList *l)
{
  return l->nelem;
}

/******************************************************************************/

void *prevElement(linkedList *l, listNode **last)
{
  // If this is the end, return null.
  if (!*last) return NULL;
  // Move the iterator along to the next element,
  // and then return the data item 
  void *e = (*last)->data;
  *last = (*last)->prev;
  return e;
}

/******************************************************************************/

void *nextElement(linkedList *l, listNode **last)
{
  // If this is the end, return null.
  if (!*last) return NULL;
  // Move the iterator along to the next element,
  // and then return the data item 
  void *e = (*last)->data;
  *last = (*last)->next;
  return e;
}

/******************************************************************************/

int linkedListContains(linkedList *l, void *e)
{
  listNode *iter = topOfLinkedList(l);
  void *this;
  
  while((this = nextElement(l, &iter)))
    if (this==e) return 1;
  
  return 0;
}

/******************************************************************************/

listNode *topOfLinkedList(linkedList *l)
{
  return l->head;
}

/******************************************************************************/
// Note: this does not free the memory associated with the removed element.

void removeFromLinkedList(linkedList *l, listNode *ln)
{
  if (!ln){
    fprintf(stderr, "Error: Tried to remove null element from linkedList.\n");
    exit(1);
  }  
  // This could be the top of the linkedList: if it is: make sure we change the
  // head pointer to the new value.
  if (ln->prev)
    ln->prev->next = ln->next;
  else
    l->head = ln->next;
 
  // This could be the bottom of the linkedList: make sure we change the last pointer.
  // if it is.   
  if (ln->next)
    ln->next->prev = ln->prev;
  else 
    l->last = ln->prev;
    
  // Free the node, and update the element count.
  push(l->deadNodes, ln);
  l->nelem --;
}

/******************************************************************************/
// This will free the elements the linkedList, along with the elemnts
// using the given destructor.
void freeLinkedList(linkedList *l, void (*destructor)(void *e))
{
  
  listNode *thisNode = l->head;
  while (thisNode)
  {
    listNode *tmp = thisNode->next;    
    if (destructor) destructor(thisNode->data);
    free(thisNode);
    thisNode = tmp;
  }
  freeStack(l->deadNodes, free);
  free(l);
}

/*******************************************************************************
* This is a simple array-based implementation of a stack, implementing
* pop, push, isEmpty, size and empty.
* We use isEmpty to tell whether we can keep popping pointers (because we could
* pop a NULL pointer to the stack). We use the empty function to reset the
* stack to zero without popping any elements. This maintains the memory
* associated with the stack.
*******************************************************************************/

stack *newStack()
{
  stack *s = malloc(sizeof(stack));
  
  s->top   = 0;
  s->slots = 2;
  s->arr = malloc(2*sizeof(void*));
  
  return s;
}

/******************************************************************************/

int stackSize(stack *s)
{
  return s->top;
}

/******************************************************************************/
// When we are storing pointers which could be null, we need to have a check 
// to see if the stack is empty.

int isEmpty(stack *s)
{
  return (s->top == 0);
}

/******************************************************************************/
// This will cause the stack to be empty: note, we have NOT freed any memory
// associated with the elements.

void emptyStack(stack *s)
{
  s->top = 0;
}

/******************************************************************************/

void  push(stack *s, void*e)
{
  if (s->top >= s->slots)
   {
      // we have to allocate more space
      s->slots *= 2;
      s->arr = realloc(s->arr, (s->slots*sizeof(void*)));
      // check that we haven't run out of memory
      if (s->arr == NULL)
      {
         fprintf(stderr, "Error: Out of Memory.\n");
         exit(1);
      }
   }
   // add the element
   s->arr[s->top] = e;
   s->top ++;
}

/******************************************************************************/

void *pop(stack *s)
{
  // If the stack is empty
  if (s->top == 0) return NULL;  
  s->top--;
  return s->arr[s->top];
}

/******************************************************************************/

void freeStack(stack *s, void (*destructor)(void *e))
{
  void *e;
  if (destructor)
    while ((e = pop(s)))
      destructor(e); 
      
  free(s->arr);
  free(s);
}

/******************************************************************************/

/* Unit Testing */
#ifdef _TEST_
#include "assert.h"

/******************************************************************************/

void testStack()
{
  stack *s = newStack(); 
  long *lt;
  long l[5] = {12e3, 982, 1201, 34e3, 3e3};

  int i;
  for (i=0;i<5;i++)
    push(s, &l[i]);
    
  while ((lt = pop(s)) && i--)
    assert(*lt = l[i]);
    
  assert(i==0);
  freeStack(s,free);
}

/******************************************************************************/

void testLinkedList()
{
  linkedList *l = newLinkedList();
  
  long l1, l2, l3;
  l1 = 1;
  l2 = 12e3;
  l3 = 2e6;
  
  assert(linkedListSize(l) == 0);
  
  addToLinkedList(l,&l1);
  addToLinkedList(l,&l2);
  addToLinkedList(l,&l3);
  
  assert(linkedListSize(l) == 3);
  long *l4 = getFromLinkedList(l,2);
  assert(*l4 == l3);
  
  // Test iteration:
  listNode *iter = topOfLinkedList(l);
  
  long *this = nextElement(l,&iter);
  assert( *this == l1 );
  this = nextElement(l,&iter);
  assert( *this == l2 );
  
  this = nextElement(l,&iter);
  assert( *this == l3 );
  
  removeFromLinkedList(l,topOfLinkedList(l));
  assert(linkedListSize(l)==2);
  
  iter = topOfLinkedList(l);
  this = nextElement(l,&iter);
  assert(*this == l2);
  removeFromLinkedList(l,iter);
  
  assert(linkedListSize(l) == 1);
  
  freeLinkedList(l, NULL);
  
}

/******************************************************************************/

int main(int argc, char **argv)
{
  // Run some test routines.
  testStack();
  testLinkedList();
  printf("Passed.\n");
}

/******************************************************************************/
#endif


