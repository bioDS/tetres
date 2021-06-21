
Summarizing trees
*****************

.. contents::

The Centroid class
==================

TODO: Change the start default to the Frechet Mean algorithm instead of last

.. code-block:: python

    class treeoclock.summary.Centroid(variation="greedy", n_cores=None,
                                      select='random', start='last')

Applying a variation of the centroid algorithm on a set of trees.


Centroid attributes
-------------------

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

    mycen.start = "random"  # Starting the algorithm from a random tree in the tree set

    cen_tree, sos = mycen.compute_centroid(TimeTreeSet())  # Computing the centroid for an empty TimeTreeSet


Centroid variations
-------------------

The variation parameter of a :class:`Centroid` has to be one in ["inc_sub", "greedy"] (TODO: Still WIP).

.. table::

   ==========================================     =========================================================================================================================
   Variation                                         Description
   ==========================================     =========================================================================================================================
    greedy                                          Computes a centroid via the greedy path and neighbourhood search. Only considering the tree with the most imporved SoS value in each iteration.
    inc_sub                                         Starts with a subsample of trees from the set, computes the greedy centroid variant and adds more trees to the subsample until all trees are part of the sample.
    WIP                                             WIP
    WIP                                             WIP
    WIP                                             WIP
   ==========================================     =========================================================================================================================




Computing the SoS
-----------------


Sampling the neighbourhood
--------------------------


Frechet Mean
============






