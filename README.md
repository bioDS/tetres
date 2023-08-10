| :warning: Package in active developement                                                                                       |
|:-------------------------------------------------------------------------------------------------------------------------------|
| If you encounter a bug or don't find documentation for something please raise an issue. |

tetres
======

This is the development version 1.0.0 of the package tetres (TimE TREe Statistics).
The submodule trees contains the data structures to work with discrete time trees, i.e. RNNI trees, and contains the distance implementation.
The submodule summary contains implementations of the centroid algorithm, a version of Sturms algorithm (Fréchet mean) and an annotaion method for RNNI trees.


Installation
------------

Download this Github repository, go into a command line and cd to the directories `.../path/tetres/{trees, judgment}/` and run the make command in each (OpenMP is required for these c extensions).

Then this package can be installed for Python with

```
pip install -e your/path/to/tetres
```

For usage check out the [documentation](https://biods.github.io/tetres/) also availabe as a [pdf](docs/build/latex/tetres.pdf).


References
----------

If you publish a paper using the summary subpacke, please cite<br>
Lars Berling, Lena Collienne, and Alex Gavryushkin<br>
**Estimating the mean in the space of ranked phylogenetic trees**<br>
*bioRxiv* 2023<br>
[https://doi.org/10.1101/2023.05.08.539790](https://doi.org/10.1101/2023.05.08.539790) 

If you publish a paper using the judgement (GR tree diagnostic) subpacke, please cite<br>
Lars Berling, Remco Bouckaert, and Alex Gavryushkin<br>
**Automated convergence diagnostic for phylogenetic MCMC analyses**<br>
*bioRxiv* 2023<br>
[N/A](nan) 

