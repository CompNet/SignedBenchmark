# File-related constants and functions.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")



OUT_FOLDER <- "out"




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.folder.path <- function(n, k, dens, prop.mispl=NA, prop.neg=NA)
{	result <- file.path(OUT_FOLDER,paste0(paste0("n=",n),paste0("_k=",k),paste0("_dens=",sprintf("%.4f",dens))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	return(result)
}