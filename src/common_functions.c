#include "common_functions.h"

#include <stdlib.h>
#include <string.h>

matrix* get_matrix(matrix* data, const char* name){
  matrix* ptr = data;
  while(ptr != NULL){
    if( (strcmp(name, ptr->name) == 0) )
      return ptr;
    ptr = ptr->next;
  }
  return NULL;
}

matrix* get_last_node(matrix* data){
  matrix* ptr = data;
  do{
    if(ptr->next == NULL)
      break;
    ptr = ptr->next;
  }while(1);

  ptr->next = malloc(sizeof(matrix));
  memset(ptr->next, 0, sizeof(matrix));
  return ptr->next;
}