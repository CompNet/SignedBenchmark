#############################################################################################
# Functions used to plot signed networks, and their node partition if any.
# 
# 11/2015 Vincent Labatut
#############################################################################################
library("igraph")





#############################################################################################
# Plot the specified signed graph, generating files of the specified formats.
#
# g: signed graph to plot.
# plot.file: base name (and path) of the generated files.
# format: format(s) to handle (among "PDF", "PNG", and NA for screen).
#
# returns: the same graph, with the spatial positions stored as nodal attributes x and y.
#############################################################################################
plot.network <- function(g, membership=NA, plot.file, format=c("PDF","PNG",NA))
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
	{	gpos <- delete.edges(graph=g,edges=which(E(g)$weight<0))
#		lyt <- layout.kamada.kawai(graph=gpos)
		lyt <- layout.fruchterman.reingold(graph=gpos)
#		lyt <- layout.circle(graph=gpos)
		# store spatial positions as nodal attributes
		V(g)$x <- lyt[,1]
		V(g)$y <- lyt[,2]
	}
	else
		lyt <- cbind(V(g)$x,V(g)$y)
	
	# process each specified format
	for(frmt in format)
	{	# set plot file name
		plot.filename <- plot.file
		if(toupper(substr(plot.filename, nchar(plot.filename)-2, nchar(plot.filename)))!=toupper(frmt))
			plot.filename <- paste0(plot.filename ,".",frmt)
		# create the file
		if(!is.na(frmt))
		{	if(frmt=="PNG")
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
