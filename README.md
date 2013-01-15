epiFit
======

This programme runs stochastic evolutionary simulations for exploring the potential impact of epistatic interactions on the fitness of populations. For backgrounds and full details please refer to:

*Hemani G, Knott S, Haley C.* **An evolutionary perspective on epistasis and the missing heritability**. *PLoS Genetics (in press)*.

### Summary

The objective here is to find how genetic variation changes over time when the causal variants are epistatic and affect fitness directly.

To this end, populations are simulated with genetic effects segregating, and we will assess:
 - How long the mutation remains polymorphic under selection
 - What the variance decomposition is (additive vs non-additive)
 - How different statistical tests perform in detecting them

We can choose to vary a number of different conditions including:
 - Population size
 - Sample size for the statistical tests
 - Which GP maps to assess
 - How many QTLs in each population
 - What is the LD between observed SNPs and causal variants
 - What are the starting allele frequencies of the causal variants
 - How many generations to run each population
 - How many repeats of each population to run

The programme is written in C, so it should run fast and will be easy to compile on Linux / Mac.

## Installation

Requires `gcc` and `gfortran`. To install (on Mac or Linux) simply clone the repo and run

    make

This will create an execultable called `epiFit`.


## How to run

`epiFit` takes just one argument, a text file that has the parameters for the simulation. An easy way to manage this is to wrap the whole thing in an R script. Please see `run_simulations.R` for step-by-step details on how to use this programme.



## Acknowledgements

My name is [Gibran Hemani][0], this work was produced as part of my PhD thesis at [The Roslin Institute][1] at the University of Edinburgh under the supervision of [Chris Haley][2] and [Sara Knott][3]

## License

epiSpaces is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

epiSpaces is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see [http://www.gnu.org/licenses/][4].

 [0]:http://www.complextraitgenomics.com/the_team/index.php#gibran_hemani
 [1]:http://www.roslin.ac.uk
 [2]:http://www.roslin.ed.ac.uk/chris-haley/
 [3]:http://www.ed.ac.uk/schools-departments/biology/evolutionary-biology/staff-profiles?id=sknott&cw_xml=homepage.php
 [4]:http://www.gnu.org/licenses/

