###############################################################################
# Main program, generates the networks and applies InfoMap (also generates a few
# plots).
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/assess.R")
source("src/generate.R")




# set up the parameters
n <- 1000									# number of nodes
k <- 5										# number of (same-sized) clusters
dens <- 0.001								# constant density
prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links




# generate the networks
generate.signed.graphs(n, k, dens, prop.mispls, prop.negs)
plot.graph.stats(n, k, dens, prop.mispls, prop.negs)




# apply InfoMap
apply.infomap(n, k, dens, prop.mispls, prop.negs)
plot.algo.stats(n, k, dens, prop.mispls, prop.negs)




###############################################################################
# TODO
###############################################################################
