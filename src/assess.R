###############################################################################
# Functions used to partition the artificially generated graphs and plot the results.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/common.R")
source("src/file.R")
source("src/plot.R")




###############################################################################
# Applies the given algorithm to the signed graphs previously generated
# using the random model. The identified partition is recorded in the same folder
# as the graph.
#
# algo.name: signed graph partitioning algo name or community detection algo name. Infomap or ExCC
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the graphs are complete).
###############################################################################
apply.partitioning.algo <- function(algo.name, n, k, dens, prop.mispl, prop.neg, network.no.list)
{	
#	if(dens==1)
#		prop.negs <- NA
	tlog(0, paste0("Start to apply ",algo.name," to the previously generalted collection of signed network",sep=""))
	tlog(2,"Parameters:")
	tlog(4,"n=",n)
	tlog(4,"k=",k)
	tlog(4,"dens=",dens)
	tlog(4,"prop.mispls=",paste(prop.mispls, collapse=", "))
	tlog(4,"prop.negs=",paste(prop.negs, collapse=", "))
	membership <- rep(1:k,each=n%/%k)
	
	for(network.no in network.no.list){
		tlog(2,"Processing network.no=",network.no)
	
		# process separately each value for the proportion of negative links
		for(prop.neg in prop.negs)
		{	tlog(2,"Processing prop.neg=",prop.neg)
			
			# process separately each value for the proportion of misplaced links
			for(prop.mispl in prop.mispls)
			{	tlog(4,"Processing prop.mispl=",prop.mispl)
				folder <- get.network.folder.path(n, k, dens, prop.mispl, prop.neg, network.no)
				
				# read the graph
				file.graph <- file.path(folder,paste0(GRAPH.FILENAME,".graphml"))
				tlog(6,"Reading the graph from ",file.graph)
				if(!file.exists(file.graph))
					tlog(2,"WARNING: Could not find graph file ",file.graph,", probably couldn't be generated due to some parameter problem.")
				else
				{	g <- read.graph(file=file.graph,format="graphml")
					
					# apply the partitioning algorithm
					tlog(6,"Applying Infomap")
					idx <- which(E(g)$weight<0)
					
					mbr = NA
					
					# ==================================================================
					# apply the partitioning algorithm
					# ==================================================================
					
					tlog(6, paste("Applying ", algo.name, sep=""))
					
					if(algo.name == IM){
						# ======== preprocessing for unsigned community detection algos
						idx <- which(E(g)$weight<0)
						g <- delete_edges(graph=g, edges=idx)
						
						res <- infomap.community(graph=g,e.weights=NULL)
						mbr <- as.vector(membership(res)) # membership function from igraph package
					} else if(algo.name == ExCC){
						# a temporary file in order to get the result of the algo
#						output.path = file.path(folder, "temp-ExCC-result.txt")
						network.path = file.path(folder, paste0(GRAPH.FILENAME,".G"))						
						mbr = run.ExCC(network.path, folder)
#						unlink(output.path) # delete the temporary file
					}
					else { # unknown algo name
						return(-1)
					}
					# =================================================================
					
					# record the result
					file.part <- file.path(folder,paste(algo.name,".txt", sep=""))
					tlog(6,"Recording the partition in file ",file.part)
					write.table(x=mbr,file=file.part,row.names=FALSE,col.names=FALSE)
					
					# plot the graph and detected partition
					file.plot <- file.path(folder, paste(algo.name,".pdf", sep=""))
					tlog(6,"Plotting the graph in file ",file.plot)
					g <- plot.network(g, membership, plot.file=file.plot, format="PDF")
				}
			}
		}
	}
	tlog(0, paste(algo.name," processing over",sep=""))
	return(0) # success
}	




###############################################################################
# Generates performance plots, to check how the partitioning algorithms are
# affected by changes in the graphs.
#
# algo.name: signed graph partitioning algo name or community detection algo name. Infomap or ExCC
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the graphs are complete).
###############################################################################
plot.algo.stats <- function(algo.name, n, k, dens, prop.mispls, prop.negs, network.no.list)
{	
#	if(dens==1)
#		prop.negs <- NA
	tlog(0,"Start to compute the algorithm statistics")
	
	CODE_IMBALANCE <- "Imbalance (proportion of misplaced links)"
	CODE_NMI <- "Normalized Mutual Information"
	CODES <- c(CODE_IMBALANCE, CODE_NMI)
	FILE_NAMES <- c()
	FILE_NAMES[CODE_IMBALANCE] <- paste0(algo.name, "-imbalance")
	FILE_NAMES[CODE_NMI] <- paste0(algo.name, "-nmi")
	gt.membership <- rep(1:k,each=n%/%k)
	
	
	res <- list()
	for(network.no in network.no.list){
	
		# load each graph and process its stats
		res[[network.no]] = list()
		
		res[[network.no]][[CODE_IMBALANCE]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		res[[network.no]][[CODE_NMI]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		for(i in 1:length(prop.mispls))
		{	for(j in 1:length(prop.negs))
			{	folder <- get.network.folder.path(n, k, dens, prop.mispls[i], prop.negs[j], network.no)
				file <- file.path(folder,paste0(GRAPH.FILENAME,".graphml"))
				if(!file.exists(file))
					tlog(2,"WARNING: Could not find graph file ",file,", probably couldn't be generated due to some parameter problem.")
				else
				{	g <- read.graph(file=file,format="graphml")
					membership <- as.matrix(read.table(file=file.path(folder,paste0(algo.name,".txt"))))	
					el <- get.edgelist(graph=g,names=FALSE)
					comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
					res[[network.no]][[CODE_IMBALANCE]][i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
					res[[network.no]][[CODE_NMI]][i,j] <- compare(gt.membership,membership,method="nmi")
				}
			}
		}
	}
	
	# generate the plots
	for(network.no in network.no.list){
		res2 = res[[network.no]] # res2 is a list containing a matrix
		
		for(code in CODES)
		{	data <- res2[[code]]
			if(!all(is.na(data)))
			{	# function of the proportion of misplaced links
				folder <- get.plot.folder.path(n, k, dens, network.no=network.no)
				dir.create(path=folder, showWarnings=FALSE, recursive=TRUE)
				plot.file <- file.path(folder, paste0(FILE_NAMES[code],"_vs_propmisp.PDF"))
				pdf(file=plot.file,bg="white")
					plot(NULL, xlim=c(min(prop.mispls),max(prop.mispls)), 
							ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)), 
							xlab="Desired proportion of misplaced links", ylab=code)
					cc <- 1
					for(j in 1:length(prop.negs))
					{	lines(x=prop.mispls, y=data[,j], col=COLORS[cc])
						cc <- cc + 1
					}
					legend(x="topright",fill=COLORS,legend=prop.negs, title="Negative links")
				dev.off()
				
				# function of the proportion of negative links
				if(dens!=1 && !all(is.na(prop.negs)))
				{	cc <- 1
					plot.file <- file.path(folder, paste0(FILE_NAMES[code],"_vs_propneg.PDF"))
					pdf(file=plot.file,bg="white")
						plot(NULL, xlim=c(min(prop.negs),max(prop.negs)), 
								ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)), 
								xlab="Desired proportion of negative links", ylab=code)
						for(i in 1:length(prop.mispls))
						{	lines(x=prop.negs, y=data[i,], col=COLORS[cc])
							cc <- cc + 1
						}
						legend(x="topright",fill=COLORS,legend=prop.mispls, title="Misplaced links")
					dev.off()
				}
			}
		}
	}
	
	tlog(0,"Computing of the algorithm statistics over")
}	




###############################################################################
# Test
###############################################################################
#n <- 1000									# number of nodes
#k <- 5										# number of clusters
#dens <- 0.001								# constant density
#prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
#prop.negs <- seq(from=0, to=1, by=0.1)		# proportion of negative links
#apply.infomap(n, k, dens, prop.mispls, prop.negs)
#plot.algo.stats(n, k, dens, prop.mispls, prop.negs)
