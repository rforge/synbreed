
# coding genotypic data
codeGeno <- function(gpData,impute=FALSE,impute.type=c("fix","random","family","beagle","beagleAfterFamily"),replace.value=NULL,maf=NULL,nmiss=NULL,label.heter="AB",keep.identical=TRUE,verbose=FALSE){

  #============================================================
  # read information from arguments
  #============================================================
  
  noHet <- is.null(label.heter) # are there only homozygous genotypes?, we need this for random imputation
  
  # check for class 'gpData'
  if(class(gpData)=="gpData"){
    if(is.null(gpData$geno)) stop("no genotypic data available")
    else data <- gpData$geno
    # family information (population structure) for genotypic data
    popStruc <- gpData$covar$family[gpData$covar$genotyped] 
    if(gpData$info$codeGeno & !noHet){
       warning("assuming heterozygous genotypes coded as 1. Use 'label.heter' to specify if that is not the case")
       label.heter <- "1"
    } 
  }
  # atm other formats are supported too
  else{         
   if(impute & impute.type %in% c("beagle","beagleAfterFamily")) stop("using Beagle is only possible for a gpData object")
   data <- gpData
   popStruc <- NULL
   gpData <- list(geno=data)
   gpData$map <- NULL
  }                                
  #  catch errors
  if(class(data)!= "data.frame" & class(data) != "matrix") stop("wrong data format")
  if (any(colMeans(is.na(data))==1)) warning("markers with no nonmissing values in data")
  
  # number of genotypes
  n <- nrow(data)

  # keep names of data object
  cnames <- colnames(data)
  rnames <- rownames(data)
  dataRaw <- matrix(unlist(data),nrow=n)

  # elements from control list
  if(impute)  impute.type <- match.arg(impute.type)
  # catch errors 
  if (impute){
  if(!is.logical(impute)) stop("impute has to be logical")
  if(impute.type=="fix" & is.null(replace.value)) stop("'replace.value' must be given for impute.type='fix'")
  # imputing with family information
  if(impute.type=="family" & is.null(popStruc)) stop(paste("family information needed, but '",substitute(gpData),"$covar$family' is empty",sep=""))
  if(impute.type=="family" & !is.null(popStruc)){
    if(length(popStruc)!=n) stop("population structure must have equal length as obsersvations in genotypic data")
    if(any(is.na(popStruc))) warning("missing values in family information, imputation is likely to be incomplete")
  }
  }
 
  #============================================================
  # step 1  - remove makrers with more than nmiss fraction of missing values (optional, argument nmiss>0)
  #============================================================
 
   if(!is.null(nmiss)){
    if(nmiss<0 | nmiss>1) stop("'nmiss' must be in [0,1]")
    which.miss <- apply(is.na(dataRaw),2,mean,na.rm=TRUE)<=nmiss 
    dataRaw <- dataRaw[,which.miss]
    if (verbose) cat("step 1 :",sum(!which.miss),"marker(s) removed with >",nmiss*100,"% missing values \n")
    cnames <- cnames[which.miss]
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.miss,]
    }
    else{
      dataRaw <- dataRaw
      if (verbose) cat("step 1 : No markers removed due to fraction of missing values \n")
    }
    
  #============================================================
  # step 2  - coding alleles
  #============================================================

   # identify heterozygous genotypes
   if(!is.null(label.heter)){
    if (is.character(label.heter)) label.heter <- label.heter # 1 label for heterozygous
    else{
      if (is.function(label.heter)){                          # multiple labels for heterozygous values
          is.heter <- label.heter
          label.heter <- unique(dataRaw[which(is.heter(dataRaw),arr.ind=TRUE)])
          
      
      }
      else stop("label.heter must be a character string or a function")
      } 
   # make sure that NA is not in label.heter
   # otherwise missing values would be masked
   label.heter <- label.heter[!is.na(label.heter)]
   }
   
   # function to recode alleles within one locus : 0 = major, 2 = minor 
   codeNumeric <- function(x){
    # names of alleles ordered by allele frequency
    alleles <-  names(table(x)[order(table(x),decreasing=TRUE)])
    # do not use heterozygous values
    alleles <- alleles[!alleles %in% label.heter]
    if (length(alleles)>2) stop("only biallelic marker allowed (check for heterozygous genotypes)")
    x[x %in% alleles] <- (as.numeric(factor(x[x %in% alleles],levels=alleles))-1)*2
    return(x)
   }

  if (verbose) cat("step 2 : Recoding alleles \n")
  # apply function on whole genotypic data
  res <- apply(dataRaw,2,codeNumeric)
 
   # set heterozygous genotypes as 1
   res[res %in% label.heter] <- 1
   res <- matrix(as.numeric(res),nrow=n)
               
  # coding of SNPs finished

  #============================================================
  # step 3  - imputing missing genotypes  (optional, argument impute=TRUE)
  #============================================================
  
  # start of imputing
  if(impute){
  
  # number of markers
   M <- ncol(res)          
  
  
   # number of missing values
   nmv <- sum(is.na(dataRaw))
 
  # initialize counter  
   cnt1 <- 0   # for nr. of imputations with family structure
   cnt2 <- 0   # for nr. of random imputations
   cnt3 <- 0   # for expected fractions of correct omputations by family structure = p^2 + (1-p)^2

  # if impute.type="fix", replace missing values according to specified value
  if(impute.type=="fix"){  
     res[is.na(res)] <- replace.value
     if (verbose) cat("step 3 : Replace missing values by",replace.value," \n")
  } 

  # impute missing values according to population structure
  if(impute.type %in% c("family" ,"beagleAfterFamily")){
   if (verbose) cat("step 3 : Imputing of missing values by family information \n")
   # initialize counter (- number of heterozygous values) 
    # loop over all markers
    for (j in 1:M){
    
        if(sum(!is.na(res[,j]))>0){
        if(j==1) ptm <- proc.time()[3]
        try({
        # compute population structure  as counts
        poptab <- table(popStruc,res[,j])
        rS <- rowSums(poptab)     
         
        # continue only if there are missing values
        if(sum(is.na(res[,j]))>0 ){
          # compute otherstatisticx
          major.allele <- unlist(attr(poptab,"dimnames")[[2]][apply(poptab,1,which.max)])
          
          # look if SNP is segregating  for this population
          polymorph <- apply(poptab,1,length) >1 & ( apply(poptab,1,min) != 0)
         
          # count missing vlalues
          nmissfam <- tapply(is.na(res[,j]),popStruc,sum)
          
          # must be a named list
          names(major.allele) <- names(polymorph)
          
          # loop over all families          
          for ( i in rownames(poptab)[nmissfam>0] ){                          
            # impute values for impute.type="family" : all missing genotypes
            if (impute.type=="family" | (impute.type=="BeagleAfterfamily" & is.na(gpData$map$pos[j]) )){
              res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(ifelse(polymorph[i],sample(c(0,2),size=nmissfam[i],prob=c(0.5,0.5),replace=TRUE),rep(major.allele[i],nmissfam[i])))
              # update counter
              ifelse(polymorph[i],cnt2 <- cnt2 + nmissfam[i],cnt1 <- cnt1 + nmissfam[i])  
            }
            if (impute.type=="beagleAfterFamily"){
              # impute values for impute.type="beagleAfterfamily"  : only monomorph markers
              res[is.na(res[,j]) & popStruc == i ,j] <- as.numeric(ifelse(polymorph[i],NA,rep(major.allele[i],nmissfam[i])))
              # update counter
              ifelse(polymorph[i],cnt2 <- cnt2 + 0,cnt1 <- cnt1 + nmissfam[i])  
            }  
          }
        }
        if(j==ceiling(M/100)) cat("         approximative run time ",(proc.time()[3] - ptm)*99," seconds ... \n",sep=" ")
     }) # end try
    }   # end of if(sum(!is.na(res[,j]))>0)
    }  # end of marker loop
    
  # run beagle
  }
  if(impute.type %in% c("beagle","beagleAfterFamily")){
   
   if (verbose) cat("step 3 : Imputing of missing values by Beagle \n")
     chr <- unique(gpData$map$chr)
     markerTEMP <- gpData
     markerTEMP$geno <- res
     rownames(markerTEMP$geno) <- rownames(gpData$geno)
     colnames(markerTEMP$geno) <- rownames(gpData$map) 
    
     # loop over chromosomses
     for (lg in seq(along=chr)){
      
      sel <- rownames(gpData$map[is.na(gpData$map$pos) | gpData$map$chr!= chr[lg],])
      if (length(sel)>0) markerTEMPbeagle <- discard.markers(markerTEMP,which=sel)
      else markerTEMPbeagle <- markerTEMP    # this occurs for only 1 chr

      # recode for Beagle
      markerTEMPbeagle$geno[markerTEMPbeagle$geno==0] <- "AA"
      markerTEMPbeagle$geno[markerTEMPbeagle$geno==1] <- "AB"
      markerTEMPbeagle$geno[markerTEMPbeagle$geno==2] <- "BB"
     
      # update counter
      cnt2 <- cnt2 + sum(is.na(markerTEMPbeagle$geno))
     
      # write input files for beagle
      pre <- paste("chr",chr[lg],sep="")
      # create new directory "beagle" for beagle input and output files
      if(!"beagle" %in% list.files()) system("mkdir beagle")
      write.beagle(markerTEMPbeagle,file.path(getwd(),"beagle"),prefix=pre)
      system("cd beagle")
      system(paste("java -Xmx1000m -jar beagle.jar unphased=beagle/",pre,"input.bgl markers=beagle/",pre,"marker.txt missing=NA out=",sep=""))
      system("cd..")
      system(paste("gzip -d -f beagle/",pre,"input.bgl.dose.gz",sep=""))
      
      # read data from beagle
      resTEMP <- read.table(paste("beagle/",pre,"input.bgl.dose",sep=""),header=TRUE,row.names=1)
      resTEMP <- t(resTEMP[,-c(1:2)])
      
      # convert dose to genotypes
      if(noHet){
        resTEMP[resTEMP<1] <- 0
        resTEMP[resTEMP>=1] <- 2
      }
      else{
         resTEMP <- round(resTEMP,0) # 0, 1, and 2
      }
      
      if (length(sel)>0) res[,!rownames(gpData$map) %in% sel] <- resTEMP
      else res <- resTEMP
    }  
   } 
    

   # impute missing values with no population structure
    if(impute.type=="random"){
    if (verbose) cat("step 3 : Random imputing of missing values \n")
       # initialize counter (- number of heterozygous values) 
        for (j in 1:M){
             cnt2 <- cnt2 + sum(is.na(res[,j]))
             # estimation of running time after the first iteration
             if(j==1) ptm <- proc.time()[3]

              p <- mean(res[,j],na.rm=TRUE)/2  # minor allele frequency
              cnt3 <- cnt3 + sum(is.na(res[,j])) * (p^2+(1-p)^2)
              if(noHet){        # assuming only 2 homozygous genotypes
                  res[is.na(res[,j]),j] <- sample(c(0,2),size=sum(is.na(res[,j])),prob=c(1-p,p),replace=TRUE)
              }
              else{                            # assuming 3 genotypes
                  res[is.na(res[,j]),j] <- sample(c(0,1,2),size=sum(is.na(res[,j])),prob=c((1-p)^2,2*p*(1-p),p^2),replace=TRUE)
              }
             if(j==ceiling(M/100) & verbose) cat("         approximate run time ",(proc.time()[3] - ptm)*99," seconds \n",sep=" ")
        }   
   }    
  
  #============================================================
  # step 4 - recoding
  #============================================================
  
  
  # recode again if allele frequeny changed to to imputing
    if(any(colMeans(res,na.rm=TRUE)>1)){
      if (verbose) cat("step 4 : Recode alleles due to imputation \n")
      res[,which(colMeans(res,na.rm=TRUE)>1)] <- 2 - res[,which(colMeans(res,na.rm=TRUE)>1)]     
    }
    else{
     if (verbose) cat("step 4 : No recoding of alleles necessary after imputation \n") 
     }
    }
  
  #============================================================
  # step 5 - remove markers with minor allele frequency < maf  (optional, argument maf>0)
  #============================================================
  
  if(!is.null(maf)){
    if(maf<0 | maf>1) stop("'maf' must be in [0,1]")
    which.maf <- apply(res,2,mean,na.rm=TRUE)>=2*maf 
    if (verbose) cat("step 5 :",sum(!which.maf),"marker(s) removed with maf <",maf,"\n")
    res <- res[,which.maf]
    cnames <- cnames[which.maf] 
    # update map
    if(!is.null(gpData$map)) gpData$map <- gpData$map[which.maf,]
  }
  else {
   if (verbose) cat("step 5 : No markers discarded due to minor allele frequency \n")
  }
  
  #============================================================
  # step 6 - discard duplicated markers   (optional, argument keep.identical=FALSE)
  #============================================================
  
  if(!keep.identical){
       which.duplicated <- duplicated(res,MARGIN=2)
       res <- res[,!which.duplicated]
       if (verbose) cat("step 6 :",sum(which.duplicated),"duplicated marker(s) removed \n")
       cnames <- cnames[!which.duplicated]
       # update map 
       if(!is.null(gpData$map)) gpData$map <- gpData$map[!which.duplicated,]
  }
  else{
  if (verbose) cat("step 6 : No duplicated markers discarded \n")
     }
     
  #============================================================
  # step 7 - restoring original data format
  #============================================================
  
   #if (verbose) cat("step 7 : Restoring original data format \n")
   if(is.matrix(data)){
     res <- matrix(res,nrow=n)
     rownames(res) <- rnames
     colnames(res) <- cnames
   }
   if(is.data.frame(data)){
      res  <- as.data.frame(res)
      colnames(res) <- cnames
      row.names(res) <- rnames
   }

  #============================================================
  # print summary of imputation
  #============================================================
  
  if(impute){
    cat("\n")
    cat("Summary of imputation \n")
    cat(paste("  total number of missing values                :",nmv,"\n"))
    if(impute.type %in% c("family","beagleAfterFamily")) cat(paste("  number of imputations by family structure     :",cnt1,"\n"))
    if(impute.type %in% c("beagle","beagleAfterFamily")) cat(paste("  number of Beagle imputations                  :",cnt2,"\n"))
    else cat(paste("  number of random imputations                  :",cnt2,"\n"))
    if(impute.type=="family") cat(paste("  approximate fraction of correct imputations   :",round((cnt1+0.5*cnt2)/(cnt1+cnt2),3),"\n"))
    if(impute.type=="random") cat(paste("  approximate fraction of correct imputations   :",round((cnt3)/cnt2,3),"\n"))
  }
  
  # overwrite original genotypic data
  if(class(gpData)=="gpData") {
    gpData$geno <- res
    gpData$info$codeGeno <- TRUE
  }
  else gpData <- res

  # return a gpData object (or a matrix)
  return(gpData)
}


