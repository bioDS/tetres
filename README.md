# README

tree.c reads a tree file, computes a shortest path with FINDPATH between the first two trees of that file, and returns all trees on this shortest path in "output/fp.rtree".

## Ranked tree file format:
1st line: number of leaves
2nd line: number of trees in file
3rd and following lines: trees in list of sets representation

### Example for an input file:

> 5
> 2
> [{3,5},{2,4},{2,3,4,5},{1,5,4,3,2}]
> [{1,2},{1,2,3},{1,2,3,4},{1,2,3,4,5}]
