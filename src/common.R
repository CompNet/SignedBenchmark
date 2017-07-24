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
