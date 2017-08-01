###############################################################################
# Main program, generates the networks and applies InfoMap (also generates a few
# plots).
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/define-consts.R")
source("src/assess.R")
source("src/generate.R")




# set up the parameters
n <- 30									# number of nodes
k <- 3										# number of (same-sized) clusters
dens <- 0.05								# constant density
prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links




# generate the networks
generate.signed.graphs(n, k, dens, prop.mispls, prop.negs)
plot.graph.stats(n, k, dens, prop.mispls, prop.negs)


# apply InfoMap
res = apply.partitioning.algo(IM, n, k, dens, prop.mispls, prop.negs)
if(res == -1) tlog("unknown partitioning algo name")
plot.algo.stats(IM, n, k, dens, prop.mispls, prop.negs)

# apply ExCC
res = apply.partitioning.algo(ExCC, n, k, dens, prop.mispls, prop.negs)
if(res == -1) tlog("unknown partitioning algo name")
plot.algo.stats(ExCC, n, k, dens, prop.mispls, prop.negs)





###############################################################################
# TODO
# - generate plots of the raw graphs, and of the detected partitions as well (use the script from netvotes)
# - ExCC is really slow. try improving the actual code source in order to run faster (Ask Zacarie for lazy constraint)
# - record execution times for both partitioning methods
###############################################################################