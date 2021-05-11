
Working with time trees
***********************

.. contents::

The TimeTree class
==================

A :class:`TimeTree` can be initialized with a given newick string using the `ete3.Tree <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#trees>`_ constructor and its `format options <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`_.
Additionally, a :class:`TREE` object (:ref:`c classes`) is generated and saved in the :class:`TimeTree` and used for efficient distance computations.

TimeTree attributes
-------------------

.. table::

   ========================================     ========================================================================================================
   Method                                       Description
   ========================================     ========================================================================================================
     :attr:`TimeTree.etree`                     returns the :class:`ete3.Tree` object
     :attr:`TimeTree.ctree`                     returns the respective :class:`TREE` object
     :attr:`len(TimeTree)`                      returns the number of leaves of the :class:`TimeTree`
     :attr:`TimeTree.fp_distance(t)`            returns the findpath distance to another :class:`TimeTree` *t*
     :attr:`TimeTree.fp_path(t)`                returns a *list* of :class:`TREE` objects
     :attr:`TimeTree.get_newick(format)`        returns the *write()* function of the :class:`ete3.Tree` with the specified *format*, defaults to *format=5*
     :attr:`TimeTree.copy()`                    returns a deep copy of the current :class:`TimeTree` (specifically used to generate a new :class:`TREE` object)
     :attr:`TimeTree.neighbours()`              returns a list of :class:`TimeTree`'s containing all neighbours at distance *1*
     :attr:`TimeTree.rank_neighbours()`         returns a list of :class:`TimeTree`'s containing only neighbours one *rank move* away
     :attr:`TimeTree.nni_neighbours()`          returns a list of :class:`TimeTree`'s containing only neighbours one *NNI move* away
   ========================================     ========================================================================================================

This is an example of how to access the different attributes of a TimeTree object:

.. code-block:: python

    from treeoclock.trees.time_trees import TimeTree


    # Initialize a time tree from a newick string
    tt = TimeTree("((1:3,5:3):1,(4:2,(3:1,2:1):1):2);")

    tt.ctree  # the TREE class object

    tt.etree  # the ete3.Tree object

    len(tt)  # Number of leaves in the tree tt --> 5

    tt.fp_distance(tt)  # Distance to another TimeTree --> 0

    tt.fp_path(tt)  # A shortest path to another TimeTree --> []

    tt.get_newick()  # Returns the newick string in ete3 format=5

    ttc = tt.copy()  # ttc contains a deep copy of the TimeTree tt

    tt.neighbours()  # a list of TimeTree objects each at distance one to tt

    tt.rank_neighbours()  # list of TimeTree obtained by doing all possible rank moves on tt

    tt.nni_neighbours()  # list of TimeTree obtained by doing all possible NNI moves on tt

ete3 functionalities
--------------------

Via the :class:`ete3.Tree` object the respective function of the :mod:`ete3` package are available for a :class:`TimeTree` object.
For example drawing and saving a tree to a file:

.. code-block:: python

    from treeoclock.trees.time_trees import TimeTree

    tt = TimeTree("((1:3,5:3):1,(4:2,(3:1,2:1):1):2);")

    # Automatically save the tree to a specific file_path location
    tt.etree.render('file_path_string')

    # Defining a layout to display internal node names in the plot
    def my_layout(node):
        if node.is_leaf():
            # If terminal node, draws its name
            name_face = ete3.AttrFace("name")
        else:
            # If internal node, draws label with smaller font size
            name_face = ete3.AttrFace("name", fsize=10)
        # Adds the name face to the image at the preferred position
        ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    ts.show_branch_length = True
    ts.show_scale = False

    # Will open a separate plot window, which also allows interactive changes and saving the image
    tt.etree.show(tree_style=ts)

See the :mod:`ete3` `documentation <http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html>`_ for more options.

The TimeTreeSet class
=====================

The TimeTreeSet class is an iterable list of TimeTree objects.


include the class documentation herer

Reading Trees
-------------

A TimeTreeSet object can be initialized with a path to a nexus file.


.. code-block:: python

    from treeoclock.trees import time_trees

    tts = TimeTreeSet(path_to_nexus_file.nex)

    for tree in tts:
        # tree is a TimeTree object
        print(tree.get_newick(f=9))
    tts[0]  # trees are accessible via the index
    len(tts)  # Returns the number of trees in the TimeTreeSet object



Writing trees
-------------
Still WIP

Random trees
============

STILL WIP


## Random tree generation via ete3

## Random tree set generation




Comparing trees
===============

Findpath
++++++++

To calculate the distance between two time trees it is possible to either just compute the distance as an integer or
to get a shortest path of time trees as a list of TREE(MISSING LINK and still WIP) objects.

.. code-block:: python

    from treeoclock.trees.findpath_distance import findpath_distance

    from treeoclock.trees.findpath_path import findpath_path






Combining multiple TimeTreeSets
===============================

Still WIP

More functions
==============

.. _c classes:

Classes for the c library
-------------------------

These are found in the _ctrees.py module.

The TREE and NODE class are the tree representation in the c code.

The TREELIST class is the output when computing a shortest path with the findpath algorithm in c.


Class converter functions
-------------------------

These are found in the _converter.py module.



The CTREE doc
=============

# Data structures

## Node
- int parent -- index of parent
- int children[2] -- index of children
- int time -- Time of the node
- Each internal node is assigned a unique time. All leaves have times of 0. If the node is a leaf node, the children are -1.

## Tree
- Node * tree -- The nodes of the tree. Length of 2*n-1.
- leaves:
    - in first n positions of array
    - children of leaves: -1, int parent = n + rank of parent - 1
- internal nodes:
    - in last n-1 positions of array => internal node of rank r has index n + r - 1
    - parent of root: -1
- double * leaves -- Leaf lengths of the tree.
- double scale -- Scaling value used to convert DCT tree back to a time tree.
- long num_leaves -- Number of leaves.
- long root_time -- Time of the root node.
- long sos_d -- Sum of squared distances of the tree to each tree in a given tree list. Used to summarize trees. Otherwise -1.

## Tree_List
- Tree * -- list of trees (as described above)
- int num_trees -- length of the tree list

## Paths
- long moves[path_length][2] -- The path matrix.
- long length -- The distance value of the path.
- Each row path[i] represents a move (i + 1)th move on path
- in first column: index of rank of lower node bounding interval of move: path[i][0] = k if move i happens on interval [k,k+1]. If the move was a length move, this column represents the node index of the move.
- in second column: path[i][1] in {0,1,2,3,4} ; 0 if rank move, 1 if NNI move where children[0] stays child of lower node, 2 if NNI move where children[1] stays child of lower node, 3 if length move where k.time was increased, 4 if length move where k.time was decreased.

## Centroid List
- char ** tree_strings -- List of trees represented as Strings.
- double ** leavesarray -- The leaf lengths of each tree.
- long * distances -- The sum of squared distances of each tree.
- double * scales -- The scaling value of each tree, used to convert trees back to time trees.
- long length -- The length of the arrays.