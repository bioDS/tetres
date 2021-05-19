
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

A :class:`TimeTreeSet` is an iterable list of :class:`TimeTree` objects, which can be initialized with a nexus file.

========================================     ========================================================================================================
   Method                                       Description
========================================     ========================================================================================================
:attr:`TimeTreeSet.map`                      a dictionary conataining the taxa to integer translation from the nexus file
:attr:`TimeTreeSet.trees`                    a list of :class:`TimeTree` objects
:attr:`TimeTreeSet[i]`                       returns the :class:`TimeTree` at :attr:`TimeTreeSet.trees[i]`
:attr:`len(TimeTreeSet)`                     returns the number of trees in the list :attr:`TimeTreeSet.trees`
:attr:`TimeTreeSet.fp_distance(i, j)`        returns the distances between the trees at postition i and j
:attr:`TimeTreeSet.fp_path(i, j)`            returns a shortest path (list of :class:`TREE`) between the trees at postition i and j
========================================     ========================================================================================================

Reading Trees
-------------

A TimeTreeSet object can be initialized with a path to a nexus file.

.. code-block:: python

    from treeoclock.trees.time_trees import TimeTreeSet


    # Initializing with a path to a nexus tree file
    tts = TimeTreeSet("path_to_nexus_file.nex")

    tts.map  # a dictionary keys:int and values:string(taxa)

    tts.trees  # A list of TimeTree objects

    for tree in tts:
        # tree is a TimeTree object
        ...
    tts[0]  # trees are accessible via the index

    len(tts)  # Returns the number of trees in the TimeTreeSet object

    tts.fp_distance(i, j)  # Returns the distance between trees i and j
    tts.fp_path(i, j)  # Returns a shortest path between trees i and j


Writing trees
-------------
Still WIP

Random trees
============

STILL WIP


Combining multiple TimeTreeSets
===============================

Still WIP

General Functions
=================

A list of the direct functions and their arguments.

=============================================    =====================================================================================
   Function                                       Description
=============================================    =====================================================================================
:attr:`time_trees.neighbourhood(tree)`              returns a list of :class:`TimeTree` objects containing the one-neighbours of tree
:attr:`time_trees.get_rank_neighbours(tree)`        returns a list of :class:`TimeTree` objects containing the rank neighbours of tree
:attr:`time_trees.get_nni_neighbours(tree)`         returns a list of :class:`TimeTree` objects containing the NNI neighbours of tree
:attr:`time_trees.read_nexus(file)`                 returns a list of :class:`TimeTree` objects contained in given the nexus file
:attr:`time_trees.get_mapping_dict(file)`           returns a dictionary containg the taxa to integer transaltion of the given file
:attr:`time_trees.findpath_distance(t1, t2)`        Computes the distance between t1 and t2
:attr:`time_trees.findpath_path(t1, t2)`            Computes the path between t1 and t2
=============================================    =====================================================================================

.. note::
    Both functions :attr:`time_trees.findpath_distance(t1, t2)` and :attr:`time_trees.findpath_path(t1, t2)`
    can be called with t1 and t2 being either a :class:`TREE`, :class:`TimeTree` or :class:`ete3.Tree`


.. _c classes:

Classes for the c library
=========================

These classes are found in the :file:`_ctrees.py` module.
The corresponding CDLL c library is generated from :file:`findpath.c`.

NODE
----

- :attr:`parent`: index of the parent node (int, defaults to -1)
- :attr:`children[2]`: index of the two children ([int], defaults to [-1, -1])
- :attr:`time`: Time of the node (int, defaults to 0)

.. note::
    The attribute :attr:`time` is currently not being used!

TREE
----

- :attr:`num_leaves`: Number of leaves in the tree (int)
- :attr:`tree`: Points to a :class:`NODE` object (POINTER(:class:`NODE`))
- :attr:`root_time`: Time of the root :class:`Node` (int)

.. note::
    The attribute :attr:`root_time` is currently not being used!

TREELIST
--------

- :attr:`num_trees`: Number of trees in the list (int)
- :attr:`trees`: List of trees (POINTER(:class:`TREE`))

Class converter functions
=========================

These are found in :file:`_converter.py`.

===========================================     ================================================================
   Function                                       Description
===========================================     ================================================================
    :attr:`_converter.ete3_to_ctree(tree)`      traverses an :class:`ete3.Tree` and construct the correct :class:`TREE`
    :attr:`_converter.ctree_to_ete3(ctree)`     recursively traverses a :class:`TREE` and generates an :class:`ete3.Tree`
===========================================     ================================================================