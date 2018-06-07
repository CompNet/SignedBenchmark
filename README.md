
SignedBenchmark v1
==================
*Benchmark to study partitioning problems on signed graphs*

* Copyright 2017 Vincent Labatut 

SignedBenchmark is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see `licence.txt`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/SignedBenchmark
* Contact: Vincent Labatut <vincent.labatut@univ-avignon.fr>

-----------------------------------------------------------------------

# Description
This set of R scripts was designed to randomly generate signed graphs possessing some form of community structure,
in order to assess partitioning algorithms. 


# Data


# Organization
Here are the folders composing the project:
* Folder `src`: contains the source code (R scripts).
* Folder `out`: contains the files produced by our scripts.


# Installation
1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`expm`](https://cran.r-project.org/web/packages/expm/index.html): required for certain signed graph layouts (tested with version 	0.999-2).
   * [`igraph`](http://igraph.org/r/): required (tested with version 1.0.1).
     * first, type in terminal for XML package: sudo apt-get install libxml2-dev
     * in R, type: install.packages("XML")
     * then, type in R: install.packages("igraph")
3. Download this project from GitHub and unzip the archive.
4. Install IBM Cplex 12.8.0
   * for ubuntu, type the following command:
     * sudo ./cplex_studio12.8.0.linux-x86-64.bin 
       * the default installation location for education version is: /opt/ibm/ILOG/CPLEX_Studio128 
       * the default installation location for trial version is:  /opt/ibm/ILOG/CPLEX_Studio_Community128/cplex/bin/x86-64_linux/
       * If you use trial version or want to change the default installation folder, update the corresponding variable "CPLEX.BIN.PATH" located in "define-consts.R". Otherwise, no need to update it
5. Install Open Java 1.8 (our ExCC jar file is compiled with that - **TODO**)


# Use
In order to replicate the experiments from the article, perform the following operations:

1. Open the `R` console.
2. Set the current projetct directory as the working directory, using `setwd("my/path/to/the/project/SignedBenchmark")`.
3. Run `src/main.R`


# Extension


# Dependencies
* [`igraph`](http://igraph.org/r/) package: used to build and handle graphs.
* [`expm`](https://cran.r-project.org/web/packages/expm/index.html) package: power of matrices.

# To-do List
* generate plots of the raw graphs, and of the detected partitions as well (use the script from netvotes)
* ExCC is too slow. Try improving the actual code source. (Ask Zacarie for lazy constraint)
* Record execution times for both partitioning methods once we improved ExCC

# References
* N/A
