
Working with time trees
***********************

.. contents::

The TimeTree class
==================

This is the main tree class in this python package.
TimeTrees can be initialized with a given newick string of the form MISSING.
This newick is then passed on to the ete3 Tree() constructor to generate an ete3 tree (MISSING LINK).
Afterwards this ete3 tree will additionally be converted to a different representation (MISSING LINK) for efficient
distance computation.

Allows to convert a tree into a newick string, either via the get_newick(f=5) function which infers the ete3 write() function
or via the write_newick(file_path, f=5) function which writes the newick string into the file file_path with the write(out_file=file_path, format=f)
function.


include the class documentation here for all the functions of a TimeTree object

.. code-block:: python

    from treeoclock.trees import time_tree

    t = TimeTree(MISSING_NWKSTRING)

    # All the funcitons for a TimeTree ?!

    len(t)  # Number of leaves in the tree t
    t.get_newick(f=9)  # Returns the newick string in ete3 format=9

ete3 functionalities
--------------------



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





Annotate a tree with a treeset
==============================

Still WIP

Combining multiple TimeTreeSets
===============================

Still WIP

More functions
==============

Classes for the c library
-------------------------

These are found in the _ctrees.py module.

The TREE and NODE class are the tree representation in the c code.

The TREELIST class is the output when computing a shortest path with the findpath algorithm in c.


Class converter functions
-------------------------

These are found in the _converter.py module.

also Something about the 'private' funcitons here