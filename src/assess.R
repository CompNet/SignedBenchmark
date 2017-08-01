###############################################################################
# Functions used to partition the artificially generated graphs and plot the results.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/define-consts.R")
source("src/common.R")
source("src/file.R")
source("src/plot.R")


########################################################################
# Extract membership vector from the partition result file
#
# file.name: the partition result file that ExCC outputs/creates
########################################################################
load.ExCC.partition <- function(file.name)
{	# open and read the file
	con <- file(file.name, "r")
	lines <- readLines(con)
	close(con)
	
	# process the file content
	i <- 1
	line <- lines[i]
	res <- list()
	
	# TODO: change here if the result file has more information than just the partition
	# in that case, put this line: while(line!="")
	while(!is.na(line) && substr(line, 1, 1) == "[") # line!=""
	{  # process current line
		line <- strsplit(x=line, "[", fixed=TRUE)[[1]][2]
		line <- strsplit(x=line, "]", fixed=TRUE)[[1]][1]
		
		# we increment by 1 at the end because C++ starts counting from 0
		nodes <- as.integer(strsplit(x=line,", ", fixed=TRUE)[[1]]) + 1
		res[[length(res)+1]] <- nodes
		
		# process next line
		i <- i + 1
		line <- lines[i]  
	}
	
	
	# build the membership vector
	mx <- max(unlist(res))
	membership <- rep(NA,mx)
	for(i in 1:length(res))
	{  nodes <- res[[i]]
		membership[nodes] <- i 
	}
	
#	# record the partition using the internal format
#	write.table(x=membership, file=partition.file, row.names=FALSE, col.names=FALSE)
	
	return(membership)
}


########################################################################
# Run ExCC (method exact based on CC problem), save the partition result into a file,
# then convert partition result into membership vector
# (each line corresponds to a node and indicates cluster id).
#
# network.path: file path of signed graph whose the extension is .G (Mario's graph input type)
# output.path: result of ExCC which basically shows detected partition
########################################################################
run.ExCC = function(network.path, output.path){
	
	cmd = 
		paste(
			"java",
#			
			paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
			"-jar",
			ExCC.JAR.PATH,
			network.path,
			">", # redirects output into a file
			output.path,
			sep=" "
		)
	
	print(cmd)
	system(cmd, wait=TRUE)	
	
	membership = load.ExCC.partition(output.path)
	return(membership)
}



# ==============================================================================
# ==============================================================================

###############################################################################
# Applies the InfoMap and ExCC algorithms to the signed graphs previously generated
# using the random model. The identified partition is recorded in the same folder
# as the graph.
#
# algo.name: signed graph partitioning algo name or community detection algo name
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
apply.partitioning.algo <- function(algo.name, n, k, dens, prop.mispls, prop.negs)
{	tlog(0, paste("Start to apply ",algo.name," to the previously generalted collection of signed network",sep=""))
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
			folder <- get.folder.path(n, k, dens, prop.mispl, prop.neg)
			
			# read the graph
			file.graph <- file.path(folder,"network.graphml")
			tlog(6,"Reading the graph from ",file.graph)
			g <- read.graph(file=file.graph,format="graphml")
			

			# ==================================================================
			# apply the partitioning algorithm
			# ==================================================================
			
			tlog(6, paste("Applying ", algo.name, sep=""))
			
			if(algo.name == IM){
				# ======== preprocessing for unsigned community detection algos
				g = remove.neg.edges.from.graph(g) # remains only positive edges
				
				res <- infomap.community(graph=g,e.weights=NULL)
				mbr <- as.vector(membership(res)) # membership function from igraph package
			} else if(algo.name == ExCC){
				# a temporary file in order to get the result of the algo
				output.path = file.path(folder, "temp-ExCC-result.txt") 
				signed.g.file.path = file.path(folder, "signed.G") 
				
				mbr = run.ExCC(signed.g.file.path, output.path)
				unlink(output.path) # delete the temporary file
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
	
	tlog(0, paste(algo.name," processing over",sep=""))
	return(0) # success
}	




###############################################################################
# Generates performance plots, to check how the partitioning algorithms are
# affected by changes in the graphs.
#
# algo.name : algo name
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network.
###############################################################################
plot.algo.stats <- function(algo.name, n, k, dens, prop.mispls, prop.negs)
{	
	CODE_IMBALANCE <- "Imbalance (proportion of misplaced links)"
	CODE_NMI <- "Normalized Mutual Information"
	CODES <- c(CODE_IMBALANCE, CODE_NMI)
	FILE_NAMES <- c()
	FILE_NAMES[CODE_IMBALANCE] <- paste0(algo.name, "-imbalance")
	FILE_NAMES[CODE_NMI] <- paste0(algo.name, "-nmi")
	gt.membership <- rep(1:k,each=n%/%k)
	
	# load each graph and process its stats
	res <- list()
	res[[CODE_IMBALANCE]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	res[[CODE_NMI]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
	for(i in 1:length(prop.mispls))
	{	for(j in 1:length(prop.negs))
		{	folder <- get.folder.path(n, k, dens, prop.mispls[i], prop.negs[j])
			g <- read.graph(file=file.path(folder,"network.graphml"),format="graphml")
			membership <- as.matrix(read.table(file=file.path(folder,paste0(algo.name,".txt"))))
			el <- get.edgelist(graph=g,names=FALSE)
			comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
			res[[CODE_IMBALANCE]][i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
			res[[CODE_NMI]][i,j] <- compare(gt.membership,membership,method="nmi")
		}
	}
	
	# generate the plots
	for(code in CODES)
	{	data <- res[[code]]
		
		# function of the proportion of misplaced links
		plot.file <- file.path(get.folder.path(n, k, dens), paste0(FILE_NAMES[code],"_vs_propmisp.PDF"))
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.mispls),max(prop.mispls)), 
				ylim=c(min(data),max(data)), 
				xlab="Desired proportion of misplaced links", ylab=code)
		cc <- 1
		for(j in 1:length(prop.negs))
		{	lines(x=prop.mispls, y=data[,j], col=COLORS[cc])
			cc <- cc + 1
		}
		legend(x="topright",fill=COLORS,legend=prop.negs, title="Negative links")
		dev.off()
		
		# function of the proportion of negative links
		cc <- 1
		plot.file <- file.path(get.folder.path(n, k, dens), paste0(FILE_NAMES[code],"_vs_propneg.PDF"))
		pdf(file=plot.file,bg="white")
		plot(NULL, xlim=c(min(prop.negs),max(prop.negs)), 
				ylim=c(min(data),max(data)), 
				xlab="Desired proportion of negative links", ylab=code)
		for(i in 1:length(prop.mispls))
		{	lines(x=prop.negs, y=data[i,], col=COLORS[cc])
			cc <- cc + 1
		}
		legend(x="topright",fill=COLORS,legend=prop.negs, title="Misplaced links")
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
#apply.infomap(n, k, dens, prop.mispls, prop.negs)
#plot.algo.stats(n, k, dens, prop.mispls, prop.negs)
