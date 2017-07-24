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
# prop.mispl: proportion of misplaced links.
# prop.neg: proportion of negative links.
#
# returns: a signed graph randomly generated according to the specified parameters.
###############################################################################
generate.signed.graph <- function(n, membership, dens, prop.mispl, prop.neg)
{	# init probabilities
	p.pos.int <- (1-prop.neg)*(1-prop.mispl)*dens
	p.neg.int <- prop.neg*prop.mispl*dens
	p.none.int <- 1 - p.pos.int - p.neg.int
	p.pos.ext <- (1-prop.neg)*prop.mispl*dens
	p.neg.ext <- prop.neg*(1-prop.mispl)*dens
	p.none.ext <- 1 - p.pos.ext - p.neg.ext
	
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
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
generate.signed.graphs <- function(n, k, dens, prop.mispls, prop.negs)
{	tlog(0,"Start to generate the collection of signed network")
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
			
			# generate the graph
			tlog(6,"Generating the graph")
			g <- generate.signed.graph(n, membership, dens, prop.mispl, prop.neg)
			
			# record the graph
			folder <- get.folder.path(n, k, dens, prop.mispl, prop.neg)
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
plot.graph.stats <- function(n, k, dens, prop.mispls, prop.negs)
{	CODE_DENSITY <- "Density"
	CODE_PROP_NEG <- "Observed proportion of negative links"
	CODE_PROP_MISP <- "Observed proportion of misplaced links"
	CODE_PROP_EXT <- "Observed proportion of external links"
	CODES <- c(CODE_DENSITY, CODE_PROP_NEG, CODE_PROP_MISP, CODE_PROP_EXT)
	FILE_NAMES <- c()
	FILE_NAMES[CODE_DENSITY] <- "density"
	FILE_NAMES[CODE_PROP_NEG] <- "prop_neg_links"
	FILE_NAMES[CODE_PROP_MISP] <- "prop_misp_links"
	FILE_NAMES[CODE_PROP_EXT] <- "prop_ext_links"
	membership <- rep(1:k,each=n%/%k)
	
	# taken from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
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
	res <- list()
	res[[CODE_DENSITY]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	res[[CODE_PROP_NEG]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	res[[CODE_PROP_MISP]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	res[[CODE_PROP_EXT]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	for(i in 1:length(prop.mispls))
	{	for(j in 1:length(prop.negs))
		{	folder <- get.folder.path(n, k, dens, prop.mispls[i], prop.negs[j])
			g <- read.graph(file=file.path(folder,"network.graphml"),format="graphml")
			res[[CODE_DENSITY]][i,j] <- graph.density(graph=g)
			res[[CODE_PROP_NEG]][i,j] <- length(which(E(g)$weight<0)) / ecount(g)
			el <- get.edgelist(graph=g,names=FALSE)
			comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
			res[[CODE_PROP_MISP]][i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
			res[[CODE_PROP_EXT]][i,j] <- length(which(!comembership)) / ecount(g)
		}
	}
	
	# generate the plots
	for(code in CODES)
	{	data <- res[[code]]
		
		# function of the proportion of misplaced links
		plot.file <- file.path(get.folder.path(n, k, dens), paste0(FILE_NAMES[code],"_vs_propmisp.PDF"))
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.mispls),max(prop.mispls)), ylim=c(0,1), xlab="Desired proportion of misplaced links", ylab=code)
		c <- 0
		for(j in 1:length(prop.negs))
		{	lines(x=prop.mispls, y=data[,j], col=colors[c])
			c <- c + 1
		}
		legend(x="topright",fill=colors,legend=prop.negs, title="Negative links")
		dev.off()
		
		# function of the proportion of negative links
		c <- 0
		plot.file <- file.path(get.folder.path(n, k, dens), paste0(FILE_NAMES[code],"_vs_propneg.PDF"))
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.negs),max(prop.negs)), ylim=c(0,1), xlab="Desired proportion of negative links", ylab=code)
		for(i in 1:length(prop.mispls))
		{	lines(x=prop.negs, y=data[i,], col=colors[c])
			c <- c + 1
		}
		legend(x="topright",fill=colors,legend=prop.negs, title="Misplaced links")
		dev.off()
	}
}	




###############################################################################
# Test
###############################################################################
#n <- 1000									# number of nodes
#k <- 5										# number of clusters
#dens <- 0.001								# constant density
#prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
#prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links
#generate.signed.graphs(n, k, dens, prop.mispls, prop.negs)
plot.graph.stats(n, k, dens, prop.mispls, prop.negs)
