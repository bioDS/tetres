# Computing nearest neighbour interchange distances between ranked phylogenetic trees - Code Repository

This repository contains the code for computing the RNNI distance using the algorithm FINDPATH from the paper [Computing nearest neighbour interchange distances between ranked phylogenetic trees](https://doi.org/10.1007/s00285-021-01567-5) by Lena Collienne and Alex Gavryushkin


## Compilation

`make`

The executable program **tree** asks for a file containing ranked trees and computes the RNNI distance between the given trees.

## Files

| File			|	Description
---			|	---
| tree.c | FINDPATH algorithm |
| example.trees | example of input file for trees in the format described below |

### Input format for tree file

> number_of_leaves  
> number_of_trees  
> tree_1  
> tree_2  
> ...

The input trees need to be given in the cluster representation.