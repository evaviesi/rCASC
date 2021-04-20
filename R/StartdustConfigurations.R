
#' @title Stardust Configurations
#' @description This function executes the method StartdustPermutation 5
#'  times, each one with a different configuration for spaceWeight
#'  parameter (0,0.25,0.5,0.75,1). Requires ggplot2 installed.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param tissuePosition, file with tissue position name with extension
#' @param nPerm, number of permutations to perform the pValue to evaluate clustering
#' @param permAtTime, number of permutations that can be computes in parallel
#' @param percent, percentage of randomly selected cells removed in each permutation
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param logTen, 1 if the count matrix is already in log10, 0 otherwise
#' @param pcaDimensions, 	0 for automatic selection of PC elbow.
#' @param seed, important value to reproduce the same results with same input
#' @param sparse, boolean for sparse matrix
#' @param format, output file format csv or txt

#' @author Giovanni Motterle, giovanni [dot] motterle [at] studenti [dot] univr [dot] it, University of Verona
#'
#' @return 5 results directories corrispondent to the 5 method configurations and a violin
#'  plots comparison of the CSSs of each configuration
#' @export

StartdustConfigurations <- function(group=c("sudo","docker"), scratch.folder, 
  file, tissuePosition, nPerm, permAtTime, percent, separator, logTen=0, 
  pcaDimensions=5, seed=1111, sparse=FALSE, format="NULL"){

  library(ggplot2)
  if(!file.exists(file))
    stop("data file don't exists")
  if(!file.exists(tissuePosition))
    stop("tissue positions file don't exists")
  if(!dir.exists(scratch.folder))
    stop("scratch folder don't exists")

  if(!sparse){
    data.folder=dirname(file)
  }
  else {
    data.folder=paste(strsplit(dirname(file),"/")[[1]][-length(strsplit(dirname(file),"/")[[1]])],collapse="/")
  }

  profileDistance=2
  spotDistance=2

  zero <- ""
  zerodue <- ""
  zerocinque <- ""
  zerosette <- ""
  uno <- ""
  for(i in c(0,0.25,0.5,0.75,1)){
    StardustPermutation(group,scratch.folder,file, tissuePosition, 
      spaceWeight = i, nPerm, permAtTime, percent, separator,
      logTen, pcaDimensions, seed, sparse, format)
    cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
    cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))
    permAnalysisSeurat(group="docker",scratch.folder = scratch.folder,file=file, nCluster=cluster,separator="\t",sp=0.8)

    current.data.folder = paste(data.folder,"/",i,sep="")
    system(paste("mkdir -p ",current.data.folder,sep=""))
    system(paste("cp -r ",data.folder,"/Results/* ",current.data.folder,sep=""))
    system(paste("rm -rf ",data.folder,"/Results",sep=""))
    path = paste(current.data.folder,"/filtered_expression_matrix/",cluster,"/filtered_expression_matrix_scoreSum.txt",sep="")
    if(i == 0){
      zero <- read.delim(path, header=FALSE, row.names=1)
    }
    else if (i == 0.25) {
      zerodue <- read.delim(path, header=FALSE, row.names=1)
    }
    else if (i == 0.5) {
      zerocinque <- read.delim(path, header=FALSE, row.names=1)
    }
    else if (i == 0.75) {
      zerosette <- read.delim(path, header=FALSE, row.names=1)
    }
    else if (i == 1) {
      uno <- read.delim(path, header=FALSE, row.names=1)
    }
  }
  system(paste("rm ",data.folder,"/*.log",sep=""))
  system(paste("rm ",data.folder,"/containers.txt",sep=""))
  system(paste("rm ",data.folder,"/ExitStatusFile",sep=""))
  system(paste("rm ",data.folder,"/run.info",sep=""))
  vectorLength = dim(zero)[1]
  x = factor(c(rep("0:1",vectorLength),rep("1:4",vectorLength),rep("1:2",vectorLength),
    rep("3:4",vectorLength),rep("1:1",vectorLength)),levels=c("0:1","1:4","1:2","3:4","1:1"))
  df = data.frame(x=x,y=c(zero$V2,zerodue$V2,zerocinque$V2,zerosette$V2,uno$V2),
                  physicalDistance = c(rep("no",vectorLength),rep("yes",vectorLength*4)))
  ggplot(df,aes(x=x, y=y, fill=physicalDistance)) + geom_violin(trim=TRUE) + 
    geom_boxplot(width=0.1, fill="white") + stat_summary(fun=mean, geom="line", aes(group=1)) + 
    stat_summary(fun=mean, geom="point")  +  labs(x="Stardust configuration - space:transcripts", y = "Cell Stability Score")
  ggsave(paste(data.folder,"/violins.pdf",sep=""), width=8, height=5)
}