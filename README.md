# GeneticGraphs
This repo contains code for generating graphs that are difficult for Laplacian linear solvers
through the use of genetic algorithms.

If you use code or results please cite:

Deweese, Kevin, and John R. Gilbert. "Evolving Difficult Graphs for Laplacian Solvers."
2018 Proceedings of the Seventh SIAM Workshop on Combinatorial Scientific Computing. Society for Industrial and Applied Mathematics, 2018.

Repo is organized  as follows:
evolution dir contains the genetic algorithm code for initial graph generation, mutation, and recombination

python dir contains some python helper scripts

solvers dir contains solver code for running linear solver on graphs and managing the evolution
Todo: Clean this code up, much of it is redundant, and too many dependencies are hardcoded.
To get an idea of how this works look at the main function in belos_kosz_evolve_ratio.cpp
At the initial generation, it reads in several matricies in a test directory, performs linear solves on them
using two different solvers, keeping only the most difficult graphs for one solver compared to another
(measured in ratio of estimated flops for solver A / estimated flops for solver B),
and removing the rest. If run at a later generation, it reads a smaller set of matrices to initialize
the next generation, and then starts performing mutation and crossover on pairs of matrices and testing
their fitness with the linear solvers, again only keeping the most difficult.

scripts dir contains bash scripts for automating large evolution trials
Todo: This folder is also a mess but contains many different attempts at automating the evolution,
generally running something like 200 generations at a time so a supervisor (that's you!) can check
the evolution progress and then continue running.

Note that for these experiments I used two external solver packages that
will need to be installed separately.
For preconditioned PCG experiments we used the Belos subpackage of Trilinos,
devloped at Sandia National Laboratories.
https://github.com/trilinos/Trilinos
For KOSZ experiments we used the implementation found in
https://github.com/serbanstan/TreePCG

