###############################################################################
# Main program, generates the networks and applies InfoMap (also generates a few
# plots).
# 
# Author: Vincent Labatut 07/2017
###############################################################################

# ----------------------------------------------------
# CONSTANTS
GRAPH.FILENAME = "signed-unweighted"

OUT.FOLDER <- "out"
NETWORKS.FOLDER <- file.path(OUT.FOLDER, "networks")
PLOTS.FOLDER <- file.path(OUT.FOLDER, "plots")
LIB.FOLDER <- "lib"
# ----------------------------------------------------

library("igraph")

source("src/define-algos.R")
source("src/assess.R")
source("src/generate.R")



# set up the parameters
#n <- 60									# number of nodes
k <- 2										# number of (same-sized) clusters
dens <- 1									# constant density
prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links (ignored if the graph is complete)

	
graph.sizes = seq(from=12, to=12, by=4)
network.no.list = seq(1, 2)

for(n in graph.sizes){

	if(dens == 1){
		# prop.negs depend only on graph size when density=1
		# because we use different model hypothesis 
		membership <- rep(1:k,each=n%/%k)
		pext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
						{	n1 <- length(which(membership==r[1]))
							n2 <- length(which(membership==r[2]))
							n1 * n2
						})) / (n*(n-1)/2)
		
		prop.negs = pext 
		# proportion of negative links = proporition of links located between clusters
	}
	
	
	# generate the networks
	generate.signed.graphs(n, k, dens, prop.mispls, prop.negs, network.no.list)
	plot.graph.stats(n, k, dens, prop.mispls, prop.negs, network.no.list)
	
	

	# apply InfoMap
	res = apply.partitioning.algo(IM, n, k, dens, prop.mispls, prop.negs, network.no.list)
	if(res == -1)
		tlog("unknown partitioning algo name:",IM)
	else
		plot.algo.stats(IM, n, k, dens, prop.mispls, prop.negs, network.no.list)
	
	
	# apply ExCC
	res = apply.partitioning.algo(ExCC, n, k, dens, prop.mispls, prop.negs, network.no.list)
	if(res == -1)
		tlog("unknown partitioning algo name:",ExCC)
	else
		plot.algo.stats(ExCC, n, k, dens, prop.mispls, prop.negs, network.no.list)
}

###############################################################################
# TODO
###############################################################################
