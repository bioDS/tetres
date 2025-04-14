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

For usage check out the [documentation](https://biods.github.io/tetres/) also availabe as a [pdf](docs/build/latex/tetres.pdf).<br>
Tutorials for specific tasks are available at [tetres-tutorials](https://github.com/Lars-B/tetres-tutorials).

References
----------

If you publish a paper using the summary subpackage, please cite<br>
Lars Berling, Lena Collienne, and Alex Gavryushkin<br>
**Estimating the mean in the space of ranked phylogenetic trees**<br>
*Bioinformatics* (2024)<br>
[https://doi.org/10.1093/bioinformatics/btae514](...)

If you publish a paper using the judgement (GR tree diagnostic) subpackage, please cite<br>
Lars Berling, Remco Bouckaert, and Alex Gavryushkin<br>
**An automated convergence diagnostic for phylogenetic MCMC analyses**<br>
*IEEE/ACM Transactions on Computational Biology and Bioinformatics* (2024)<br>
[https://doi.org/10.1109/TCBB.2024.3457875](https://doi.org/10.1101/2023.08.10.552869)<br>
