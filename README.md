# AINV-honeybees
Programs to generate pedigree files and compute the inverse of honeybee specific relationship matrices

This repository contains R scripts to generate inverses of honeybee specific relationship matrices following Brascamp & Bijma, 2014's methodology (10.3389/finsc.2023.1135187) written by Pim Brascamp.
Each folder, corresponding each to different versions, contains a manual and two main scripts: one creates a pedigree file with a suitable sequence of individuals and some other inputs for the second one, mainly creating the AINV.giv (inverse relationship matrix) file. This file can then typically be used for genetic analyses.
The R program is not well-suited for large datasets because of necessary memory and computing time.
