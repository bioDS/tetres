#ifndef TREE_H_
#define TREE_H_


typedef struct Node {
  long parent;
  long children[2];
  long time;
} Node;

typedef struct Tree{
  Node * tree;
  long num_leaves;
  long root_time;
  long sos_d;
} Tree;

typedef struct Tree_List{
  Tree * trees;
  int num_trees;
} Tree_List;


typedef struct Path{
  long ** moves;
  long length;
} Path;


#endif
