# Computing nearest neighbour interchange distances between ranked phylogenetic trees - Code Repository

This repository contains the code for computing the RNNI distance using the algorithm FINDPATH from the paper [Computing nearest neighbour interchange distances between ranked phylogenetic trees](https://doi.org/10.1007/s00285-021-01567-5) by Lena Collienne and Alex Gavryushkin


## Compilation

`make`

## Reading Trees

The submodule dct_parser allows reading trees from nexus files or as newick strings

## Functions

The following C functions (tree.c) are wrapped in Python (tree_structs.py):

| Functions			|	Description
---			|	---
| findpath_distance | returns distance for two trees, using FINDPATH |
| findpath_path | returns a shortest path for two trees, using FINDPATH |

