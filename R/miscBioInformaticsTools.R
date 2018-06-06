selectQuals<-function(left,right,...){
  #last base often low qual so ignore
  q1<-sapply(qualToInts(substring(left$qual,1,nchar(left$qual)-1)),function(xx)sum(10^(-xx/10)))
  q2<-sapply(qualToInts(substring(right$qual,1,nchar(right$qual)-1)),function(xx)sum(10^(-xx/10)))
  #less than 1 expected error in both reads
  
  # higher threshold for q2?
  #selector<-q1<1 & q2<1 & !grepl('[^ACTG]',left$seq)&!grepl('[^ACTG]',right$seq)
  selector<-q1<1 & q2<2 & !grepl('[^ACTG]',left$seq)&!grepl('[^ACTG]',right$seq)
  
  return(selector)
}

castAsMatrix <- function(FastqQuality, n) {
  groupIndexArray <- cut(seq(1:length(FastqQuality)),n,labels=FALSE)
  tempMatrix <- as(FastqQuality[groupIndexArray==1],"matrix")
  for(i in 2:n) {
    tempMatrix <- rbind(tempMatrix,
                        as(FastqQuality[groupIndexArray==i],"matrix"))
  }
  return(tempMatrix)
}

# Convert QC string to Phred scores
qcToPhred <- function(qcString) {
  sapply(strsplit(qcString,split="")[[1]],
         function(X) {
           return(as.integer(charToRaw(X))-33)
         })
}

# Expected Errors From Phred
eefp <- function(phred) {
  sum(sapply(phred,function(X) {
    10^(X/-10)
  }))
}

# Expected Errors from QC String
eefqc <- function(qcString) {
  return(sum(sapply(strsplit(qcString,split="")[[1]],
                    function(X) {
                      return(10^((as.integer(charToRaw(X))-33)/-10))
                    })))
}








