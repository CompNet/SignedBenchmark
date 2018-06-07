###############################################################################
# Main program, generates the networks and applies InfoMap (also generates a few
# plots).
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/assess.R")
source("src/generate.R")

# CONSTANTS
GRAPH.FILENAME = "signed-unweighted"

# set up the parameters
#n <- 60						# number of nodes
k <- 2						# number of (same-sized) clusters
dens <- 1					# constant density
prop.mispls <- seq(from=0, to=1, by=0.05)	# proportion of misplaced links
prop.negs <- seq(from=0, to=1, by=0.05)		# proportion of negative links (ignored if the graph is complete)
	
graph.sizes = seq(from=44, to=80, by=4)
network.no.list = seq(1, 10)

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
	

	
# plot stats'ta nasil yapmis bak
	
	# generate the networks
	generate.signed.graphs(n, k, dens, prop.mispls, prop.negs, network.no.list)
	
	# TODO update based on 'network.no.list'
#	plot.graph.stats(n, k, dens, prop.mispls, prop.negs)
	
	
	# apply InfoMap
	#apply.infomap(n, k, dens, prop.mispls, prop.negs)
	#plot.algo.stats(n, k, dens, prop.mispls, prop.negs)
}

###############################################################################
# TODO
###############################################################################
