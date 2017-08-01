#############################################################################################
# Defines common functions for all scripts.
# 
# 05/2016 Vincent Labatut
#############################################################################################




#############################################################################################
# Logs the specified message on screen, adding current date and time, and possible some
# offset (to represent the hierarchy of function calls).
#
# offset: number of "." used to represent the hierarchical level of the message.
# ...: parameters fetched to the cat function.
#############################################################################################
tlog <- function(offset=NA, ...)
{	prefix <- paste("[",format(Sys.time(),"%a %d %b %Y %X"),"] ",sep="")
	if(!is.na(offset))
	{	if(is.numeric(offset))
		{	os <- paste(rep(".",offset), sep="", collapse="")
			prefix <- paste(prefix, os, sep="")
		}
		else
			prefix <- paste(prefix, offset, sep="")
	}
	cat(prefix, ..., "\n", sep="")
}



#############################################################################################
# Removes negative edges from the graph in order to run community detection algos
#
# g: an igraph graph object
#############################################################################################
remove.neg.edges.from.graph = function(g){
	idx <- which(E(g)$weight<0)
	g <- delete_edges(graph=g, edges=idx)
	return(g)
}


##############################################################################
# convert "graphml" graph file into ".G" graph file.
# ".G" graph file is used by Mario

# An example:
#  4	5
#  1	2	-1
#  1	3	1
#  2	3	-0.6
#  2	4	0.42
#  3	4	0.1
#
###############################################################################
convert.graphml.into.mario.graph = function(g, mario.graph.filename){
	first.line = paste(vcount(g), ecount(g),sep="\t")
	write.graph(g, mario.graph.filename, format="ncol")
	lines = readLines(con = mario.graph.filename)
	lines2=sapply(lines, function(x) gsub(" ", "\t", x, fixed = TRUE))
	write.table(first.line, mario.graph.filename,row.names=F, col.names=F, quote=F)
	write.table(lines2, mario.graph.filename,row.names=F, col.names=F, append=T, quote=F)
}