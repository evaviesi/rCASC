
#' @title spaGCN Permutation
#' @description This function executes a ubuntu docker that produces a specific number of permutation to evaluate clustering 
#'  using spaGCN
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param h5matrix.name, a character string indicating the full path + name with extension (.h5) of the expression matrix
#' @param spotpositions.name, a character string indicating the full path + name with extension (.csv) of the spot positions
#' @param image.name, a character string indicating the full path + name with extension (.tif) of the tissue image
#' @param use_histology, a boolean indicating if spaGCN should use the hisological information or not
#' @param lResolution, double for the resolution of Louvain algorithm
#' @param pcaDimensions, number of PCA dimension to keep
#' @param p, a real number between 0 and representing the percentage of total expression contributed by neighborhoods
#' @param nPerm, number of permutations to perform the pValue to evaluate clustering
#' @param permAtTime, number of permutations that can be computes in parallel
#' @param percent, percentage of randomly selected cells removed in each permutation
#' @param seed, important value to reproduce the same results with same input

#' @author Giovanni Motterle, giovanni [dot] motterle [at] studenti [dot] univr [dot] it, University of Verona
#'
#' @return To write
#' @export

spaGCNPermutation <- function(group=c("sudo","docker"), scratch.folder,
  h5matrix.name, spotpositions.name, image.name, use_histology=TRUE, lResolution=0.8, 
  pcaDimensions=10,p=0.5,nPerm=80, permAtTime=8, percent=10, seed=1111){
  
  data.folder = normalizePath(dirname(h5matrix.name))
  matrixName = strsplit(h5matrix.name,"/",fixed = TRUE)[[1]]
  matrixName = matrixName[length(matrixName)]
  matrixName = strsplit(matrixName,".",fixed = TRUE)[[1]][1]

  lResolution = as.double(lResolution)
  if(lResolution < 0 || lResolution > 5){
    stop("Error: lResolution is not a double in [0,5]")
  }

  p = as.double(p)
  if(p > 1 || p < 0){
    stop("Parameter 'p' has not value between 0 and 1")
  }
  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    system("echo 2 > ExitStatusFile 2>&1")
    return(2)
  }

  #storing the position of the home folder
  home <- getwd()
  setwd(data.folder)
  #initialize status
  system("echo 0 > ExitStatusFile 2>&1")

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile 2>&1")
    setwd(home)
    return(10)
  }



  #check  if scratch folder exist
  if (!file.exists(scratch.folder)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
    system("echo 3 > ExitStatusFile 2>&1")
    setwd(data.folder)
    return(3)
  }
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=file.path(scratch.folder, tmp.folder)
  writeLines(scrat_tmp.folder,paste(data.folder,"/tempFolderID", sep=""))
  cat("\ncreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files


  dir.create(paste(scrat_tmp.folder,"/",matrixName,sep=""))
  dir.create(paste(data.folder,"/Results",sep=""))

  system(paste("cp ",h5matrix.name," ",scrat_tmp.folder,"/",sep=""))
  system(paste("cp ",spotpositions.name," ",scrat_tmp.folder,"/",sep=""))
  system(paste("cp ",image.name," ",scrat_tmp.folder,"/",sep=""))

  h5matrixfile = strsplit(h5matrix.name,"/",fixed = TRUE)[[1]]
  h5matrixfile = h5matrixfile[length(h5matrixfile)]
  spotpositionsfile = strsplit(spotpositions.name,"/",fixed = TRUE)[[1]]
  spotpositionsfile = spotpositionsfile[length(spotpositionsfile)]
  imagefile = strsplit(image.name,"/",fixed = TRUE)[[1]]
  imagefile = imagefile[length(imagefile)]
  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/dockerID -v ",scrat_tmp.folder,
    ":/scratch -v ", data.folder, 
    ":/data -d docker.io/giovannics/spagcnpermutation Rscript /home/main.R ",
    h5matrixfile," ",spotpositionsfile," ",imagefile," ",use_histology," ",lResolution," "
    ,pcaDimensions," ",p," ",nPerm," ",permAtTime," ",percent," ",seed," ",sep="")

  resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(data.folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("elapsed run time mins ",ptm[3]/60, sep="")

    writeLines(tmp.run,"run.info")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/dockerID", sep=""), warn = FALSE)
  system(paste("docker logs ", substr(container.id,1,12), " &> ",data.folder,"/", 
    substr(container.id,1,12),".log", sep=""))
  system(paste("docker rm ", container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  system(paste("cp -r ",scrat_tmp.folder,"/* ",data.folder,"/Results",sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  system(paste("rm -R ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),
    "containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
} 