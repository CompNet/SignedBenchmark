###############################################################################
# Functions used to partition the artificially generated graphs and plot the results.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/common.R")
source("src/file.R")




###############################################################################
# Applies the Infomap algorithm to the signed graphs previously generated
# using the random model. The identified partition is recorded in the same folder
# as the graph.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
apply.infomap <- function(n, membership, dens, prop.mispl, prop.neg)
{	tlog(0,"Start to apply infomap to the previously generalted collection of signed network")
	tlog(2,"Parameters:")
	tlog(4,"n=",n)
	tlog(4,"k=",k)
	tlog(4,"dens=",dens)
	tlog(4,"prop.mispls=",paste(prop.mispls, collapse=", "))
	tlog(4,"prop.negs=",paste(prop.negs, collapse=", "))
	membership <- rep(1:k,each=n%/%k)
	
	# process separately each value for the proportion of negative links
	for(prop.neg in prop.negs)
	{	tlog(2,"Processing prop.neg=",prop.neg)
		
		# process separately each value for the proportion of misplaced links
		for(prop.mispl in prop.mispls)
		{	tlog(4,"Processing prop.mispl=",prop.mispl)
			
			# read the graph
			tlog(6,"Reading the graph")
			folder <- get.folder.path(n, k, dens, prop.mispl, prop.neg)
			g <- read.graph(file=file.path(folder,"network.graphml"),format="graphml")
			
			# apply the partitioning algorithm
			tlog(6,"Applying Infomap")
			idx <- which(E(g)$weight<0)
			g <- delete_edges(graph=g, edges=idx)
			res <- infomap.community(graph=g,e.weights=NULL)
			mbr <- as.vector(membership(res))
			
			# record the result
			write.table(x=mbr,file=file.path(folder,"infomap.txt"),row.names=FALSE,col.names=FALSE)
		}
	}
	
	tlog(0,"Generation over")
}	




###############################################################################
# Generates performance plots, to check how the partitioning algorithms are
# affected by changes in the graphs.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
plot.algo.stats <- function(n, k, dens, prop.mispls, prop.negs)
{	# taken from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
	colors <- c(
		rgb(228,26,28,maxColorValue=255),
		rgb(55,126,184,maxColorValue=255),
		rgb(77,175,74,maxColorValue=255),
		rgb(152,78,163,maxColorValue=255),
		rgb(255,127,0,maxColorValue=255),
		rgb(255,255,51,maxColorValue=255),
		rgb(166,86,40,maxColorValue=255),
		rgb(247,129,191,maxColorValue=255),
		rgb(153,153,153,maxColorValue=255),
		rgb(0,0,0,maxColorValue=255)
	)
	
	# load each graph and process its stats
	imbalance <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	for(i in 1:length(prop.mispls))
	{	for(j in 1:length(prop.negs))
		{	folder <- get.folder.path(n, k, dens, prop.mispls[i], prop.negs[j])
			g <- read.graph(file=file.path(folder,"network.graphml"),format="graphml")
			membership <- as.matrix(read.table(file=file.path(folder,"infomap.txt")))
			el <- get.edgelist(graph=g,names=FALSE)
			comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
			imbalance[i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
		}
	}
	
	# generate the plots
	{	# function of the proportion of misplaced links
		plot.file <- file.path(get.folder.path(n, k, dens), "infomap-imbalance_vs_propmisp.PDF")
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.mispls),max(prop.mispls)), 
#				ylim=c(0,1),
				ylim=c(min(imbalance),max(imbalance)),
				xlab="Desired proportion of misplaced links", ylab="Imbalance (proportion of misplaced links)")
		c <- 0
		for(j in 1:length(prop.negs))
		{	lines(x=prop.mispls, y=imbalance[,j], col=colors[c])
			c <- c + 1
		}
		legend(x="topright",fill=colors,legend=prop.negs, title="Negative links")
		dev.off()
		
		# function of the proportion of negative links
		c <- 0
		plot.file <- file.path(get.folder.path(n, k, dens), "infomap-imbalance_vs_propneg.PDF")
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.negs),max(prop.negs)), 
#				ylim=c(0,1), 
				ylim=c(min(imbalance),max(imbalance)),
				xlab="Desired proportion of negative links", ylab="Imbalance (proportion of misplaced links)")
		for(i in 1:length(prop.mispls))
		{	lines(x=prop.negs, y=imbalance[i,], col=colors[c])
			c <- c + 1
		}
		legend(x="topright",fill=colors,legend=prop.negs, title="Misplaced links")
		dev.off()
	}
}	




###############################################################################
# Test
###############################################################################
n <- 1000									# number of nodes
k <- 5										# number of clusters
dens <- 0.001								# constant density
prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links
#apply.infomap(n, k, dens, prop.mispls, prop.negs)
plot.algo.stats(n, k, dens, prop.mispls, prop.negs)

# TODO
# - generate plots of the raw graphs, and of the detected partitions as well (use the script from netvotes)
# - pb: the proportions of pos/neg and well/mispl do not respect the parameters
