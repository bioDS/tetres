
Summarizing trees
*****************

.. contents::

The Centroid class
==================

.. code-block:: python

    class treeoclock.summary.Centroid(variation="greedy_omp", n_cores=None, select='random', start='FM', subsample_size=200,
                 tree_log_file="", max_iterations=None)


This is used to setup a Centroid object which then takes a :class:`TimeTreeSet` as input to compute the centroid summary tree.


Variation
---------

The variation parameter of a :class:`Centroid` has to be one in ["inc_sub", "greedy"] (TODO: Still WIP).

.. table::

   ==========================================     =========================================================================================================================
   Variation                                         Description
   ==========================================     =========================================================================================================================
    :ref:`var greedy`                               Computes a centroid via the greedy path and neighbourhood search. Only considering the tree with the most imporved SoS value in each iteration.
    :ref:`var incsub`                               Starts with a subsample of trees from the set, computes the greedy_omp centroid variant and adds more trees to the subsample until all trees are part of the sample.
    :ref:`var itersub`                              Starts with a subsample of trees from the set, computes the greedy_omp centroid variant and then resamples a new subset, using the previous centroid as the starting tree.
    :ref:`var separate`                             Only computes rank move neighbours if the tree contains all common clades of the tree set
    :ref:`var onlyone`                              Prefers either NNI or Rank moves and switches this if a local optimum is reached
    :ref:`var update-with-one`                      Similar to the incsub variation, only one tree at a time is added to the subsample
    :ref:`var online`                               Mimicks an online approach where samples arrive one after another and the centroid is computed after each sample starting from the previous centroid
   ==========================================     =========================================================================================================================

.. _var greedy:

Greedy
``````

.. code-block:: python

    treeoclock.summary.Centroid(variation="greedy_omp")  # default, using multiple processes in c!
    treeoclock.summary.Centroid(variation="greedy")  # pure python version


.. _var incsub:

Inc_sub
```````
The parameter subsample_size defines the size of the subsample of trees that is added each iteration.
The parameter max_iterations defines the number of iterations, if it is None the regular break is defined whenever an
iteration is not successful at improving the previous centroid.
If it is an integer then it defines the number of iterations that will subsample, if it is 0 the start tree will be returned.

.. code-block:: python

    treeoclock.summary.Centroid(variation="inc_sub", subsample_size=500, max_iterations=None)


.. _var itersub:

Iter_sub
`````````

The parameter subsample_size defines the size of the subsample of trees that is sampled each iteration.
The parameter max_iterations defines the number of iterations, if it is None the regular break is defined whenever an
iteration is not successful at improving the previous centroid.
If it is an integer then it defines the number of iterations that will subsample, if it is 0 the start tree will be returned.


.. code-block:: python

    treeoclock.summary.Centroid(variation="iter_sub", subsample_size=500, max_iterations=None)


.. _var separate:

Separate
`````````

Will only use one move, current implementation is for NNI moves only, needs to be switched in source code (_variations.py, line 147).

.. code-block:: python

    treeoclock.summary.Centroid(variation="separate")



.. _var onlyone:

Onlyone
```````

Will always do one move (starting with rank moves as of current implementation) and switch the move type whenever a local optimum is found.

.. code-block:: python

    treeoclock.summary.Centroid(variation="onlyone")


.. _var update-with-one:

update-with-one
```````````````

Similar to the inc-sub variation but only one new tree is added in each iteration.

.. code-block:: python

    treeoclock.summary.Centroid(variation="update_with_one")


.. _var online:

Online
``````

Mimicks an online approach where the trees arrive one by one in the given order.

.. code-block:: python

    treeoclock.summary.Centroid(variation="online")



Selecting a tree
----------------

This is only the case if multiple trees have the same SoS value.
The defualt is random and the options are either random, first or last.
The second two options are depending on the ordering which is dictated by the way the neighbourhood
is computed.

Starting tree
-------------

There are the options to start with the last, the first or any given index of tree from the given tree set.
The default option however is the sorted Frechet Mean tree (ref), see the doc on FM for more detail.


Subsample size
--------------

This is used by some variations and can be set to any integer number (default is 200).
This number indicates the size of the subsample that the variation will use in its iterations.
See the incsub or itersub variations

Maximal iterations
------------------

This is used to limit the number of iterations the iterative subsampling and increasing subsampling centroid versions are computing.
If it is None (default) then the regular break points of those variations apply, otherwise it will only compute upto max_iteration
many iterations.


Computing the SoS
-----------------

The n_cores parameters defines the number of cores to use, if -1 all available cores are used (default).


Tree logfile
------------

This option will write the trees of each centroid iteration to the given file path.
This includes the actual centroid as the last tree.
Can be used for further analysis.

Note that for incsub for example the tree is logged after an iteration on the subsample.
This results in much smaller log files.


Annotation of a centroid
========================

To keep the discrete ranks of a centroid use this annotation method.
Each rank get assigned the average height of that rank in the given tree set, guaranteed to keep the same ranked tree after the annotation.

.. code-block:: python

    from treeoclock.summary.annotate_centroid import annotate_centroid
    cen, sos = Centroid(variation="greedy").compute_centroid(my_tts)  # centroid of the TimeTreeSet my_tts
    annotated_cen = annotate_centroid(cen, my_tts)  # Annotation with the branch lengths from the TimeTreeSet my_tts
    # the annotated_cen is a TimeTree object for further use such as writing the newick to a file


Frechet Mean
============






