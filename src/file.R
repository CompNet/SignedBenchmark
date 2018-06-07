# File-related constants and functions.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no: network no with the same parameter sets
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.network.folder.path <- function(n, k, dens, prop.mispl=NA, prop.neg=NA, network.no=NA)
{	result <- file.path(NETWORKS.FOLDER,paste0(paste0("n=",n),paste0("_k=",k),paste0("_dens=",sprintf("%.4f",dens))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	if(!is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
	
	return(result)
}


###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no: network no with the same parameter sets
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.plot.folder.path <- function(n, k, dens, prop.mispl=NA, prop.neg=NA, network.no=NA)
{	result <- file.path(PLOTS.FOLDER,paste0(paste0("n=",n),paste0("_k=",k),paste0("_dens=",sprintf("%.4f",dens))))
	if(!is.na(prop.mispl))
		result <- file.path(result, paste0("propMispl=",sprintf("%.4f",prop.mispl)))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	if(!is.na(network.no))
		result <- file.path(result, paste0("network=", network.no))
	
	return(result)
}
