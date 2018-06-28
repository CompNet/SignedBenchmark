###############################################################################
# Functions used to randomly generate the signed graphs. They are based on
# a generalization of the simple model described in:
# 		Community structure in social and biological networks 
# 		M. Girvan & M. E. J. Newman
# 		Proceedings of the National Academy of Sciences, 2002, 99(12):7821-7826
#		DOI:10.1073/pnas.122653799
#  
# See the article for a description of how we generalized the model to
# signed networks.
# 
# Author: Vincent Labatut 07/2017
###############################################################################
library("igraph")

source("src/common.R")
source("src/file.R")
source("src/plot.R")




###############################################################################
# Builds a complete signed graph according to the specified parameters.
#
# membership: integer vector, each value represents the cluster containing the corresponding node.
# prop.mispl: proportion of misplaced links (among existing links).
#
# returns: a signed graph randomly generated according to the specified parameters.
###############################################################################
generate.complete.signed.graph <- function(membership, prop.mispl)
{	n <- length(membership)
	m <- n*(n-1)/2
	
	# generate links and weights
	weights <- rep(NA, m)
	el <- matrix(NA, nrow=m, ncol=2)
	r <- 1
	for(i in 1:(n-1))
	{	for(j in (i+1):n)
		{	el[r,] <- c(i,j)
			if(membership[i]==membership[j])
				weights[r] <- 1
			else
				weights[r] <- -1
			r <- r + 1
		}
	}
	
	# check parameters
	half.misp <- floor(m*prop.mispl/2)
	idx.p <- which(weights>0)
	idx.n <- which(weights<0)
	half.upper.bound <- min(length(idx.p), length(idx.n)) / m
	upper.bound = half.upper.bound * 2
	if(prop.mispl>upper.bound)
		stop("Parameter prop.mispl (",sprintf("%.4f",prop.mispl),") must be smaller or equal to ",sprintf("%.4f",upper.bound))
	
	# set the misplaced links
	idx.misp.p <- sample(x=idx.p, size=half.misp, replace=FALSE)
	idx.misp.n <- sample(x=idx.n, size=half.misp, replace=FALSE)
	idx.misp <- c(idx.misp.p, idx.misp.n)
	weights[idx.misp] <- -weights[idx.misp]
	
	# build graph
	g <- make_empty_graph(n,directed=FALSE)
	g <- add_edges(graph=g, edges=c(t(el)), attr=list(weight=weights))
	
	return(g)
}




###############################################################################
# Builds an incomplete signed graph according to the specified parameters.
#
# membership: integer vector, each value represents the cluster containing the corresponding node.
# dens: density of the graph.
# prop.mispl: proportion of misplaced links (among existing links).
# prop.neg: proportion of negative links (among existing links).
#
# returns: a signed graph randomly generated according to the specified parameters.
###############################################################################
generate.incomplete.signed.graph <- function(membership, dens, prop.mispl, prop.neg)
{	# init proportions
	n <- length(membership)
	qm <- prop.mispl
	qw <- 1 - qm
	tlog(8,"n=",n," qm=",qm," qw=",qw)
	# init sign probas
	p.neg <- prop.neg * dens
	p.pos <- dens - p.neg
	tlog(8,"p.neg=",p.neg," p.pos=",p.pos," (total=",p.neg+p.pos,")")
	# init position probas
	p.int <- sum(sapply(1:max(membership), function(c)
		{	n <- length(which(membership==c))
			n*(n-1)/2
		})) / (n*(n-1)/2)
	p.ext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
		{	n1 <- length(which(membership==r[1]))
			n2 <- length(which(membership==r[2]))
			n1 * n2
		})) / (n*(n-1)/2)
	tlog(8,"p.int=",p.int," p.ext=",p.ext," (total=",p.int+p.ext,")")
	# check prop.neg interval
	if(prop.mispl==0.5)
	{	tlog(8,"prop.neg (",sprintf("%.4f",prop.neg),")")
		upper.bound <- 2 * min(p.ext, p.int)
		tlog(8,"density upper bound bound: ",sprintf("%.4f",upper.bound))
		if(dens>upper.bound)
			stop("Parameter dens (",sprintf("%.4f",dens),") must be smaller or equal to ",sprintf("%.4f",upper.bound))
	}
	else
	{	if(prop.mispl<0.5)
		{	lower.bound <- max(0, (dens*(1-qm) - (1 - p.ext)) / (dens*(1-2*qm)))
			upper.bound <- min(1, (p.ext - dens*qm) / (dens*(1-2*qm)))
		}
		else if(prop.mispl>0.5)
		{	lower.bound <- max(0, (p.ext - dens*qm) / (dens*(1-2*qm)))
			upper.bound <- min(1, (dens*(1-qm) - (1 - p.ext)) / (dens*(1-2*qm)))
		}
		tlog(8,"prop.neg (",sprintf("%.4f",prop.neg),") bounds: [",sprintf("%.4f",lower.bound)," ; ",sprintf("%.4f",upper.bound),"]")
		if(prop.neg<lower.bound | prop.neg>upper.bound)
			stop("Parameter prop.neg (",sprintf("%.4f",prop.neg),") must be in [",sprintf("%.4f",lower.bound)," ; ",sprintf("%.4f",upper.bound),"]")
	}
	# init internal probas
	p.neg.int <- qm * p.neg / p.int
	p.pos.int <- qw * p.pos / p.int
	p.none.int <- 1 - p.pos.int - p.neg.int
	tlog(8,"Internal probas: neg=",sprintf("%.4f",p.neg.int)," pos=",sprintf("%.4f",p.pos.int)," none=",sprintf("%.4f",p.none.int))
	# init external probas
	p.neg.ext <- qw * p.neg / p.ext
	p.pos.ext <- qm * p.pos / p.ext
	p.none.ext <- 1 - p.pos.ext - p.neg.ext
	tlog(8,"External probas: neg=",sprintf("%.4f",p.neg.ext)," pos=",sprintf("%.4f",p.pos.ext)," none=",sprintf("%.4f",p.none.ext))
	
	# draw links
	el <- t(combn(x=n, m=2, simplify=TRUE))
	comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
	norm.int <- length(which(comembership))
	norm.ext <- length(which(!comembership))
	tlog(8,"Maximal link numbers: total=",nrow(el)," internal=",norm.int," external=",norm.ext)
	tlog(8,"Expected link numbers: total=",round(nrow(el)*dens)," positive=",round(nrow(el)*dens*(1-prop.neg))," negative=",round(nrow(el)*dens*prop.neg)," well-placed=",round(nrow(el)*dens*(1-prop.mispl))," misplaced=",round(nrow(el)*dens*prop.mispl))
	weights <- rep(NA,nrow(el))
	weights[comembership] <- sample(x=c(-1,+1,0),size=length(which(comembership)),replace=TRUE,prob=c(p.neg.int,p.pos.int,p.none.int))
	obs.p.neg.int <- length(which(weights[comembership]<0))/norm.int
	obs.p.pos.int <- length(which(weights[comembership]>0))/norm.int
	obs.p.none.int <- length(which(weights[comembership]==0))/norm.int
	tlog(8,"Drawn internal probas: neg=",sprintf("%.4f",obs.p.neg.int)," pos=",sprintf("%.4f",obs.p.pos.int)," none=",sprintf("%.4f",obs.p.none.int))
	weights[!comembership] <- sample(x=c(-1,+1,0),size=length(which(!comembership)),replace=TRUE,prob=c(p.neg.ext,p.pos.ext,p.none.ext))
	obs.p.neg.ext <- length(which(weights[!comembership]<0))/norm.ext
	obs.p.pos.ext <- length(which(weights[!comembership]>0))/norm.ext
	obs.p.none.ext <- length(which(weights[!comembership]==0))/norm.ext
	tlog(8,"Drawn external probas: neg=",sprintf("%.4f",obs.p.neg.ext)," pos=",sprintf("%.4f",obs.p.pos.ext)," none=",sprintf("%.4f",obs.p.none.ext))
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
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the density is 1).
###############################################################################
generate.signed.graphs <- function(n, k, dens=1, prop.mispls, prop.negs=NA, network.no.list)
{	
#	if(dens==1)
#		prop.negs <- NA
		
	tlog(0,"Start to generate the collection of signed networks")
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
			
			for(network.no in network.no.list){
				
				# generate the graph
				tlog(6,"Generating the graph")
				if(dens<1)
					g <- tryCatch(
						generate.incomplete.signed.graph(membership, dens, prop.mispl, prop.neg),
						error=function(e) NA
					)
				else
					g <- tryCatch(
							generate.complete.signed.graph(membership, prop.mispl),
							error=function(e) NA
					)
				
				if(!all(is.na(g)))
				{	# possibly create the output folder
					folder <- get.network.folder.path(n, k, dens, prop.mispl, prop.neg, network.no)
					dir.create(folder,showWarnings=FALSE,recursive=TRUE)
					
					# plot the graph
					file.plot <- file.path(folder,GRAPH.FILENAME)
					tlog(6,"Plotting the graph in file ",file.plot)
					g <- plot.network(g, membership, plot.file=file.plot, format="PDF", method="fruchterman.reingold")
					
					# record the graph
					file.graph <- file.path(folder,paste0(GRAPH.FILENAME,".graphml"))
					tlog(6,"Recording the graph in file ",file.graph)
					write_graph(graph=g,file=file.graph,format="graphml")
					
					# ========================================================
					# NEW
					
					# export using a format compatible with pILS
					t <- get.edgelist(graph=g) - 1	# start numbering nodes at zero
					t <- cbind(t,E(g)$weight)		# add weights as the third column
					file.graph <- file.path(folder,paste0(GRAPH.FILENAME,".G"))
					write.table(data.frame(vcount(g),nrow(t)), file=file.graph, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) # write header
					write.table(t, file=file.graph, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE) # write proper graph
				
				
					# export using a format compatible with JENA
#					n = vcount(g)
					mtrx = as_adj(graph=g,attr="weight")
					adj.mtrx.text = ""
					for(i in 1:(n-2)){
						for(j in 1:n){
							if(i<j){
								if(i == (j-1)){
									adj.mtrx.text = paste0(adj.mtrx.text, mtrx[i,j])
								} else{
									adj.mtrx.text = paste(adj.mtrx.text, mtrx[i,j])
								}
							} 
							
						}
						adj.mtrx.text = paste(adj.mtrx.text,"\n",sep="")
					}
					adj.mtrx.text = paste(adj.mtrx.text, mtrx[n-1,n], sep="")
					
					text=paste(as.character(n),"\n",sep="")
					for(v in 1:n){
						text=paste(text,(v-1),"\n",sep="")
					}
					write(x=paste0(text,adj.mtrx.text),file=file.path(folder,paste0(GRAPH.FILENAME,".jena")))
					
					# export using a format compatible with PAJEK
					write.graph(graph=g, file=file.path(folder,paste0(GRAPH.FILENAME,".net")), format="pajek")
					
					# NEW - END
					# ========================================================
				}
				else
				{	qm <- prop.mispl
					p.ext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
					{	n1 <- length(which(membership==r[1]))
						n2 <- length(which(membership==r[2]))
						n1 * n2
					})) / (n*(n-1)/2)
					
					if(qm==0.5)
					{	upper.bound <- 2 * min(p.ext, 1 - p.ext)
						tlog(8,"Could not generate the graph: the dens parameter (",sprintf("%.4f",dens),") is greater than ",sprintf("%.4f",upper.bound))
					}
					else
					{	if(qm<0.5)
						{	lower.bound <- max(0, (dens*(1-qm) - (1 - p.ext)) / (dens*(1-2*qm)))
							upper.bound <- min(1, (p.ext - dens*qm) / (dens*(1-2*qm)))
						}
						else
						{	lower.bound <- max(0, (p.ext - dens*qm) / (dens*(1-2*qm)))
							upper.bound <- min(1, (dens*(1-qm) - (1 - p.ext)) / (dens*(1-2*qm)))
						}
						tlog(6,"Could not generate the graph: the prop.neg parameter (",sprintf("%.4f",prop.neg),") is out of range [",sprintf("%.4f",lower.bound)," ; ",sprintf("%.4f",upper.bound),"]")
					}
				}
			}
			
		}
	}
	
	tlog(0,"Generation over")
}




###############################################################################
# Generates graph plots, to check if their structure corresponds to our expectations.
#
# n: number of nodes in the graph.
# k: number of clusters in the graph.
# dens: total density of the graph (counting both negative and positive links).
# prop.mispls: vector of proportions of misplaced links.
# prop.negs: vector of proportions of negative links in the network (ignored if
#			 the density is 1, i.e. we want to generate complete graphs).
###############################################################################
plot.graph.stats <- function(n, k, dens=1, prop.mispls, prop.negs, network.no.list)
{	if(dens==1)
#		prop.negs <- NA
	tlog(0,"Start to compute the graph statistics")
	
	CODE_DENSITY <- "Density"
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
	
	res <- list()
	# load each graph and process its stats
	for(network.no in network.no.list){
		res[[network.no]] = list()
		
		res[[network.no]][[CODE_DENSITY]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		res[[network.no]][[CODE_PROP_NEG]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		res[[network.no]][[CODE_PROP_MISP]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		res[[network.no]][[CODE_PROP_EXT]] <- matrix(NA,nrow=length(prop.mispls),ncol=length(prop.negs))
		
		for(i in 1:length(prop.mispls))
		{	for(j in 1:length(prop.negs))
			{	
				folder <- get.network.folder.path(n, k, dens, prop.mispls[i], prop.negs[j], network.no)
				file <- file.path(folder,paste0(GRAPH.FILENAME,".graphml"))
				if(!file.exists(file))
					tlog(2,"WARNING: Could not find graph file ",file,", probably couldn't be generated due to some parameter problem.")
				else
				{	g <- read.graph(file=file,format="graphml")
					res[[network.no]][[CODE_DENSITY]][i,j] <- graph.density(graph=g)
					res[[network.no]][[CODE_PROP_NEG]][i,j] <- length(which(E(g)$weight<0)) / ecount(g)
					el <- get.edgelist(graph=g,names=FALSE)
					comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
					res[[network.no]][[CODE_PROP_MISP]][i,j] <- length(which((!comembership & E(g)$weight>0) | (comembership & E(g)$weight<0))) / ecount(g)
					res[[network.no]][[CODE_PROP_EXT]][i,j] <- length(which(!comembership)) / ecount(g)
				}
			}
		}		
	}
	
	
	# generate the plots
	for(network.no in network.no.list){
		res2 = res[[network.no]] # res2 is a list containing a matrix
		
		for(code in CODES)
		{	data <- res2[[code]] # data is a matrix
			if(!all(is.na(data)))
			{	
				# function of the proportion of misplaced links
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
			else
				tlog(2,"WARNING: nothing to plot (no data)")
		}
	}

	tlog(0,"Computation of the graph statistics over")
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
#plot.graph.stats(n, k, dens, prop.mispls, prop.negs)

#n <- 25
#k <- 5
#membership <- rep(1:k,each=n%/%k)
#dens <- 0.1
#prop.mispl <- 0.1
#prop.neg <- 0.3
#g <- generate.incomplete.signed.graph(membership, dens, prop.mispl, prop.neg)
#plot.network(g, membership, plot.file="sqdffsd", format=NA, method="bn_zheng")

#n <- 30
#k <- 3
#membership <- rep(1:k,each=n%/%k)
#dens <- 0.2
#prop.mispl <- 0.1
#prop.neg <- 0.3
#g <- generate.incomplete.signed.graph(membership, dens, prop.mispl, prop.neg)
#plot.network(g, membership, plot.file="sqdffsd", format=NA, method="bn_zheng")

#n <- 30
#k <- 2
#membership <- rep(1:k,each=n%/%k)
#dens <- 1
#prop.mispl <- 0.1
#g <- generate.complete.signed.graph(membership, prop.mispl)
#plot.network(g, membership, plot.file="sqdffsd", format=NA, method="fruchterman.reingold")

#n <- 100									# number of nodes
#k <- 3										# number of clusters
#dens <- 1									# constant density
#prop.mispls <- seq(from=0, to=1, by=0.1)	# proportion of misplaced links
#generate.signed.graphs(n, k, dens, prop.mispls, prop.negs=NA)
#plot.graph.stats(n, k, dens, prop.mispls, prop.negs=NA)
