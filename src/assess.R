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
# Applies the InfoMap algorithm to the signed graphs previously generated
# using the random model. The identified partition is recorded in the same folder
# as the graph.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the graphs are complete).
###############################################################################
apply.infomap <- function(n, k, dens, prop.mispl, prop.neg, network.no.list)
{	
#	if(dens==1)
#		prop.negs <- NA
	tlog(0,"Start to apply InfoMap to the previously generated collection of signed networks")
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
					g <- delete_edges(graph=g, edges=idx)
					res <- infomap.community(graph=g,e.weights=NULL)
					mbr <- as.vector(membership(res))
					
					# record the result
					file.part <- file.path(folder,"infomap.txt")
					tlog(6,"Recording the partition in file ",file.part)
					write.table(x=mbr,file=file.part,row.names=FALSE,col.names=FALSE)
					
					# plot the graph and detected partition
					file.plot <- file.path(folder,"infomap.pdf")
					tlog(6,"Plotting the graph in file ",file.plot)
					g <- plot.network(g, membership, plot.file=file.plot, format="PDF")
				}
			}
		}
	}
	tlog(0,"InfoMap processing over")
}	




###############################################################################
# Generates performance plots, to check how the partitioning algorithms are
# affected by changes in the graphs.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the graphs are complete).
###############################################################################
plot.algo.stats <- function(n, k, dens, prop.mispls, prop.negs, network.no.list)
{	
#	if(dens==1)
#		prop.negs <- NA
	tlog(0,"Start to compute the algorithm statistics")
	
	CODE_IMBALANCE <- "Imbalance (proportion of misplaced links)"
	CODE_NMI <- "Normalized Mutual Information"
	CODES <- c(CODE_IMBALANCE, CODE_NMI)
	FILE_NAMES <- c()
	FILE_NAMES[CODE_IMBALANCE] <- "infomap-imbalance"
	FILE_NAMES[CODE_NMI] <- "infomap-nmi"
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
					im.membership <- as.matrix(read.table(file=file.path(folder,"infomap.txt")))
					el <- get.edgelist(graph=g,names=FALSE)
					comembership <- sapply(1:nrow(el),function(r) im.membership[el[r,1]]==im.membership[el[r,2]])
					res[[network.no]][[CODE_IMBALANCE]][i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
					res[[network.no]][[CODE_NMI]][i,j] <- compare(gt.membership,im.membership,method="nmi")
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
