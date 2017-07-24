# File-related constants and functions.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")



out.folder <- "out"




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: density of the graph.
# p.int: probability to get an internal link.
# prop.neg: proportion of negative links in the network.
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.folder.path <- function(n, k, dens, p.int, prop.neg)
{	result <- file.path(out.folder,paste0(paste0("n=",n),paste0("_k=",k),paste0("_dens=",sprintf("%.4f",dens))),
			paste0("pInt=",sprintf("%.4f",p.int)),
			paste0("propNeg=",sprintf("%.4f",prop.neg))
	)
	return(result)
}
