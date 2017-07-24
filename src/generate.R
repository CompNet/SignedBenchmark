###############################################################################
# Functions used to randomly generate the signed graphs. They are based on
# a generalization of the simple model described in:
#
# M. Girvan & M. E. J. Newman
# Community structure in social and biological networks 
# Proceedings of the National Academy of Sciences, 2002, 99(12):7821-7826
# DOI:10.1073/pnas.122653799
#  
# See the article for a description of how we generalized the model to
# signed networks.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/common.R")
source("src/file.R")




###############################################################################
# Builds a signed graph according to the specified parameters.
#
# n: number of nodes in the graph.
# membership: integer vector, each value represents the cluster containing the corresponding node.
# dens: density of the graph.
# p.int: probability to get an internal link.
# prop.neg: proportion of negative links in the network.
#
# returns: a signed graph randomly generated according to the specified parameters.
###############################################################################
generate.signed.graph <- function(n, membership, dens, p.int, prop.neg)
{	# init probabilities
	dens.pos <- dens * (1-prop.neg)
	dens.neg <- dens * prop.neg
	p.pos.ext <- dens.pos - p.pos.int
	prop <- prop.neg/(1-prop.neg)
	p.neg.int <- p.pos.ext * prop
	p.neg.ext <- p.pos.int * prop
	p.none.int <- 1 - p.neg.int - p.pos.int
	p.none.ext <- 1 - p.neg.ext - p.pos.ext
	
	tlog(8,"Internal probas: neg=",p.neg.int," pos=",p.pos.int," none=",p.none.int)
	tlog(8,"External probas: neg=",p.neg.ext," pos=",p.pos.ext," none=",p.none.ext)
	
	# draw links
	el <- t(combn(x=n, m=2, simplify=TRUE))
	comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
	tlog(8,"Maximal link numbers: total=",nrow(el)," internal=",length(which(comembership))," external=",length(which(!comembership)))
	weights <- rep(NA,nrow(el))
	weights[comembership] <- sample(x=c(-1,+1,0),size=length(which(comembership)),replace=TRUE,prob=c(p.neg.int,p.pos.int,p.none.int))
	weights[!comembership] <- sample(x=c(-1,+1,0),size=length(which(!comembership)),replace=TRUE,prob=c(p.neg.ext,p.pos.ext,p.none.ext))
	tlog(8,"Drawn link numbers: total=",length(which(weights!=0))," positive=",length(which(weights>0))," negative=",length(which(weights<0)))
	idx <- which(weights==0)
	if(length(idx)>0)
	{	el <- el[-idx,,drop=FALSE]
		weights <- weights[-idx,drop=FALSE]
	}
	
	# build graph
	g <- make_empty_graph(n,directed=FALSE)
	g <- add_edges(graph=g, edges=c(t(el)), attr=list(weight=weights))
	
	return(g)
}




###############################################################################
# Generates a collection of signed graphs corresponding to the specified parameters.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# p.ints: vector of probabilities to get an internal link.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
generate.signed.graphs <- function(n, k, dens, p.ints, prop.negs)
{	tlog(0,"Start to generate the collection of signed network")
	tlog(2,"Parameters:")
	tlog(4,"n=",n)
	tlog(4,"k=",k)
	tlog(4,"dens=",dens)
	tlog(4,"p.ints=",paste(p.ints,collapse=", "))
	tlog(4,"prop.negs=",paste(prop.negs,collapse=", "))
	membership <- rep(1:k,each=n%/%k)
	
	# process separately each value for the proportion of negative links
	for(prop.neg in prop.negs)
	{	tlog(2,"Processing prop.neg=",prop.neg)
		
		# process separately each value for the probability to get a positive internal link
		for(p.int in p.ints)
		{	tlog(4,"Processing p.int=",p.int)
			
			# generate the graph
			tlog(6,"Generating the graph")
			g <- generate.signed.graph(n, membership, dens, p.int, prop.neg)
			
			# record the graph
			folder <- get.folder.path(n, k, dens, p.int, prop.neg)
			tlog(6,"Recording the graph in folder ",folder)
			dir.create(folder,showWarnings=FALSE,recursive=TRUE)
			write_graph(graph=g,file=file.path(folder,"network.graphml"),format="graphml")
		}
	}
	
	tlog(0,"Generation over")
}




###############################################################################
# Generates graph plots, to check if their structure corresponds to our expectations.
#
# 
###############################################################################
plot.graph.stats <- function(n, k, dens, p.ints, prop.negs)
{	# load each graph and process its stats
	
	
	# generate the plots
}	




###############################################################################
# Test
###############################################################################
n <- 1000
k <- 5
dens <- 0.001
p.ints <- seq(from=dens, to=0, by=-0.0001)
prop.negs <- seq(from=0, to=1, by=0.1)
generate.signed.graphs(n, k, dens, p.ints, prop.negs)
