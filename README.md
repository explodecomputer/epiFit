epiFit
======

This programme runs stochastic evolutionary simulations for exploring the potential impact of epistatic interactions on the fitness of populations. For backgrounds and full details please refer to:

*Hemani G, Knott S, Haley C.* **An evolutionary perspective on epistasis and the missing heritability**. *PLoS Genetics (in press)*.

## Installation

Requires `gcc` and `gfortran`. To install (on Mac or Linux) simply clone the repo and run

    make

This will create an execultable called `epiFit`.


## How to run

`epiFit` takes just one argument, a text file that has the parameters for the simulation. An easy way to manage this is to wrap the whole thing in an R script. 
