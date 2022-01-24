
Summarizing trees
*****************

.. contents::

The Centroid class
==================

.. code-block:: python

    class treeoclock.summary.Centroid(variation="greedy", n_cores=None,
                                      select='random', start='last')


This class is used to setup the Centroid and then run the algorithm on a set of trees.


Centroid variations
-------------------

The variation parameter of a :class:`Centroid` has to be one in ["inc_sub", "greedy"] (TODO: Still WIP).

.. table::

   ==========================================     =========================================================================================================================
   Variation                                         Description
   ==========================================     =========================================================================================================================
    :ref:`var greedy`                               Computes a centroid via the greedy path and neighbourhood search. Only considering the tree with the most imporved SoS value in each iteration.
    :ref:`var incsub`                               Starts with a subsample of trees from the set, computes the greedy centroid variant and adds more trees to the subsample until all trees are part of the sample.
    :ref:`var itersub`                              Starts with a subsample of trees from the set, computes the greedy centroid variant and then resamples a new subset, using the previous centroid as the starting tree.
    WIP                                             WIP
    WIP                                             WIP
   ==========================================     =========================================================================================================================

.. _var greedy:

Greedy
``````

.. code-block:: python

    treeoclock.summary.Centroid(variation="greedy")


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

The n_cores parameter is used for the SOS computation

Mention the different SoS computation options
Multithreading
Multiprocessing
Benefits and downsides of both, which is used


Tree logfile
------------

This option will write the trees of each centroid iteration to the given file path.
This includes the actual centroid as the last tree.
Can be used for further analysis.

Note that for incsub for example the tree is logged after an iteration on the subsample.
This results in much smaller log files as the "greedy part" of the iteration is not logged at all because it can lead to
a tree which is not used in the end (Not necessarily an improvement)

Centroid attributes
-------------------

TODO not sure if this is necessary at all


.. table::

   ==========================================     =========================================================================================================================
   Method                                         Description
   ==========================================     =========================================================================================================================
     :attr:`Centroid.variation`                    Accesses the variation parameter which has to be in the list of specified variations.
     :attr:`Centroid.n_cores`                      Accesses the number of Threads that will be used to compute the SoS value, using ThreadPool.
     :attr:`Centroid.select`                       Accesses the selection parameter, in case of multiple options choose a selected tree i.e. the first, last or a random tree.
     :attr:`Centroid.start`                        Accesses the start parameter, specifying at which point to start the centroid computation.
     :attr:`Centroid.compute_centroid(trees)`      Computes the centroid for a given set of trees, using the options that were selected via the other attributes of the :class:`Centroid`.
   ==========================================     =========================================================================================================================


Examples of how to set and use the attributes of a :class:`Centroid`:

.. code-block:: python

    from treeoclock.trees.time_trees import TimeTreeSet
    from treeoclock.summary.centroid import Centroid


    # Initializing an empty Centroid class
    mycen = Centroid()

    mycen.variation = "inc_sub"  # Changing the variation

    mycen.n_cores = 4  # Setting the number of Threads to use

    mycen.select = "first"  # Selection parameter set to the first tree found

    mycen.start = 6  # Starting the algorithm from the sixth tree in the set (Error if non-existant)

    # Computing the centroid for an empty TimeTreeSet
    cen_tree, sos = mycen.compute_centroid(TimeTreeSet())



Frechet Mean
============






