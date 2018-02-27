#############################################################################################
# Functions used to plot signed networks, including their node partition if any.
# 
# 11/2015 Vincent Labatut
#############################################################################################
library("igraph")
library("expm")




#############################################################################################
# Colors used in the plots.
# Taken from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9
#############################################################################################
COLORS <- c(
		rgb(228,26,28,maxColorValue=255),
		rgb(55,126,184,maxColorValue=255),
		rgb(77,175,74,maxColorValue=255),
		rgb(152,78,163,maxColorValue=255),
		rgb(255,127,0,maxColorValue=255),
		rgb(255,255,51,maxColorValue=255),
		rgb(166,86,40,maxColorValue=255),
		rgb(247,129,191,maxColorValue=255),
		rgb(153,153,153,maxColorValue=255),
		rgb(0,0,0,maxColorValue=255),
		rgb(50,50,50,maxColorValue=255)
)




#############################################################################################
# This layout for signed graphs is based on some generalization of the the Laplacian for 
# signed graph. Several ones were proposed in the following articles:
#	- Spectral Analysis of Signed Graphs for Clustering, Prediction and Visualization
#		Kunegis, J.; Schmidt, S.; Lommatzsch, A.; Lerner, J.; de Luca, E. & Albayrak, S.
#		SIAM International Conference on Data Mining, 2010, p559-570.
#	- Spectral embedding of signed networks
#		Zheng, Q. & Skillicorn, D. B.
#		SIAM International Conference on Data Mining, 2015, p55-63.
#
# g: graph to process, it must possess a weight edge parameter containing the signed weights
#    of the links. It is supposed to be undirected: if it is directed, we first colapse the 
#	 links before processing it.
# method: which generalization of the signed Laplacian to use:
#			- "kunegis": Signed Laplacian (Kunegis et al. 2010)
#			- "rw_kunegis": Random walk normalized signed Laplacian (Kunegis et al. 2010)
#			- "sn_kunegis": Symmetric normalized signed Laplacian (Kunegis et al. 2010)
#			- "sn_zheng": Simple normalized signed Laplacian (Zheng & Skillicorn 2015) 
#			- "bn_zheng": Balanced normalized signed Laplacian (Zheng & Skillicorn 2015) 
#############################################################################################
layout.signed.laplacian <- function(g, method="bn_zheng")
{	# possibly make the graph undirected
	if(is.directed(g))
	{	warning("The graph should be undirected: we collapse its links before going on")
		g <- as.undirected(graph=g, mode="collapse", edge.attr.comb="mean")
	}
	
	# weighted adjacency matrix
	W <- - as.matrix(get.adjacency(graph=g, type="both", attr="weight")) # TODO don't know why we need to take the opposite of the matrix. there's a sign mistake somewhere, but couldn't find it.
	# positive weighted adjacency matrix
	Wp <- W
	Wp[which(Wp<0)] <- 0
	# negative weighted adjacency matrix
	Wm <- -W
	Wm[which(Wm<0)] <- 0
	
	# positive strength matrix
	Sp <- matrix(0,nrow(W),ncol(W))
	diag(Sp) <- apply(Wp, 1, sum)
	# negative strength matrix
	Sm <- matrix(0,nrow(W),ncol(W))
	diag(Sm) <- apply(Wm, 1, sum)
	# total strength matrix
	S <- Sp + Sm
	# its inverse matrix (simple since it's a diagonal matrix)
	Si <- S
	diag(Si) <- 1/diag(Si)
	
	# process the laplacian
	if(method=="kunegis")
		# signed Laplacian
		L <- S - W
	else if(method=="rw_kunegis")
		# random walk normalized signed Laplacian
		L <- diag(ncol(W)) - Si %*% W
	else if(method=="sn_kunegis")
		# symmetric normalized signed laplacian
		L <- diag(ncol(W)) - S %^% (-1/2) %*% W %*% S %^% (-1/2)
	else if(method=="sn_zheng")
		# simple normalized signed Laplacian
		L <- Si %*% (Sp - W)
	else if(method=="bn_zheng")
		# balanced normalized signed Laplacian
		L <- Si %*% (Sp - Sm - W)
	else
		stop("Could not recognize the specified method '",method,"'")
	
	# use zeros for isolates (division by zero on the diagonal of S)
	L[which(is.nan(L))] <- 0
	
	# take its first two Eigenvectors
	ev <- eigen(L, symmetric=TRUE)$vec
	lay <- ev[,1:2]
#	lay <- t(ev[1:2,])
#	print(lay)
	
	return(lay)
}



#############################################################################################
# Plot the specified signed graph, generating files of the specified formats.
#
# g: signed graph to plot.
# plot.file: base name (and path) of the generated files.
# format: format(s) to handle (among "PDF", "PNG", and NA for screen).
# method: layout used (signed graphs).
#
# returns: the same graph, with the spatial positions stored as nodal attributes x and y.
#############################################################################################
plot.network <- function(g, membership=NA, plot.file, format=c("PDF","PNG",NA), method="bn_zheng")
{	# setup node parameters
	vertex.sizes <- 3
	vertex.label <- NA
	
	# setup link parameters
	if(ecount(g)>0)
	{	edge.colors <- rep(NA,ecount(g))
		edge.colors[E(g)$weight<0] <- adjustcolor("RED", alpha.f=0.5)
		edge.colors[E(g)$weight>0] <- adjustcolor("GREEN", alpha.f=0.5)
		edge.widths <- abs(E(g)$weight)*1
	}
	
	# set up node colors
	if(all(is.na(membership)))
		vertex.colors <- "SkyBlue2" # default igraph color
	else
	{	palette <- rainbow(length(unique(membership)))
		vertex.colors <- palette[membership]
	}
	
	# setup layout (only if the graph doesn't already have one)
	v.att <- list.vertex.attributes(graph=g)
	if(!("x" %in% v.att) | !("y" %in% v.att)) 
	{	# layouts for unsigned graphs
		gpos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
##		lyt <- layout.kamada.kawai(graph=gpos)
#		lyt <- layout.fruchterman.reingold(graph=gpos)
##		lyt <- layout.circle(graph=gpos)
		# layouts for signed graphs
		lyt <- layout.signed.laplacian(g, method)
		# store spatial positions as nodal attributes
		V(g)$x <- lyt[,1]
		V(g)$y <- lyt[,2]
	}
	else
		lyt <- cbind(V(g)$x,V(g)$y)
	
	# process each specified format
	for(frmt in format)
	{	# create the file
		if(!is.na(frmt))
		{	# set plot file name
			plot.filename <- plot.file
			if(toupper(substr(plot.filename, nchar(plot.filename)-2, nchar(plot.filename)))!=toupper(frmt))
				plot.filename <- paste0(plot.filename ,".",frmt)
			# handle format
			if(frmt=="PNG")
			{	png(filename=plot.filename,width=800,height=800,units="px",pointsize=20,bg="white")
			}
			else if(frmt=="PDF")
			{	pdf(file=plot.filename,bg="white")
			}
		}
		
		# create the plot
		if(ecount(g)>0)
		{	plot(g, layout=lyt, #main=g$name,
				vertex.size=vertex.sizes, vertex.label=vertex.label, vertex.color=vertex.colors,
				edge.color=edge.colors, edge.width=edge.widths
#				asp=1
			)
		}
		else
		{	plot(g, layout=lyt,
				vertex.size=vertex.sizes, vertex.label=vertex.label, vertex.color=vertex.colors,
			)
		}
		title(g$name, cex.main=0.5)
		
		# finalize plot file
		if(!is.na(frmt))
			dev.off()
	}
	
	return(g)
}




#############################################################################################
### Test using the article toy graph
#g2 <-graph.empty(6,directed=FALSE)
#g2 <- add_edges(g2,edges=c(1,2,2,4,2,6,3,4,4,6,5,6))
#E(g2)$weight <- c(1,-1,-1,1,-1,1)
#plot.network(g2, membership=membership, plot.file="sqdq", format=NA, method="bn_zheng")
