
# ===============================================================
# Exact Approach: ExCC
# ===============================================================
ExCC <- "ExCC"
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,"cplex-partition.jar",sep="/") # gaia cluster - CERI
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/" # LIA server
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/" # in my computer

# ===============================================================
# Infomap
# ===============================================================
IM = "Infomap"






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
# This version of ExCC does not output all optimal solutions
#
# network.path: file path of signed graph whose the extension is .G (Mario's graph input type)
# output.folder: result folder
########################################################################
run.ExCC = function(network.path, out.folder){
	
#	cmd = 
#		paste(
#			"java",
##			
#			paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
#			"-jar",
#			ExCC.JAR.PATH,
#			network.path,
#			">", # redirects output into a file
#			output.path,
#			sep=" "
#		)


	# An example:
	# java -Djava.library.path=/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/
	# -DinFile=data/test.G -DoutDir=. -Dcp=false -DenumAll=false -Dtilim=-1 -DlazyInBB=false
	# -DusercutInBB=false -jar exe/cplex-partition.jar
	
	cmd = 
		paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", network.path),
				paste0("-DoutDir=", out.folder),
				"-Dcp=false",
				"-DenumAll=false",
				"-Dtilim=-1",
				"-DlazyInBB=false",
				"-DuserCutInBB=false",
				"-jar",
				ExCC.JAR.PATH,
				sep=" "
		)
	
	print(cmd)
	system(cmd, wait=TRUE)	
	
	# do not change the filename. This is a requirement by ExCC
	out.filename = file.path(out.folder, "ExCC-result.txt")
	membership = load.ExCC.partition(out.filename)
	return(membership)
}
