#hyplotype as a vector of elements eg. bp or loci
#among hyplotypes, homologous chromosomes put in the same order
#if there's not any separator between chromosomes within vector, "l.chromosomes" argument needed
#l.chromosomes - vector of chromosomes' lengths (chromosomes in haplotypes' order)
#lambda - mean value of Poisson distribution of number of chiasms in each chromosome, default = 2
#chiasm.before = if TRUE chiasms are formed before sampled locus or bp - chiasm in position at beginning of each chromosome cannot exist 
#meiotic - if TRUE, returns 4 gametes like in meiosis, if FALSE, returns 1 gamete
crossing.over <- function(haplotype1,haplotype2,sep="",l.chromosomes,lambda = 2,chiasm.before = T,meiotic=F){
  if(length(haplotype1) != length(haplotype2)){
    stop("Lengths of haplotypes are not equal")
  }
  # getting boundaries between chromosomes in haplotypes
  if(sep != ""){
    if(haplotype1[length(haplotype1)]==sep & haplotype2[length(haplotype2)]==sep){
      haplotype1 <- haplotype1[1:(length(haplotype1)-1)]
      haplotype2 <- haplotype2[1:(length(haplotype2)-1)]
    }
    boundaries <- (which(haplotype1 == sep))-1
    boundaries2 <- (which(haplotype2 == sep))-1
    if((length(boundaries)==length(boundaries2))==F){
      stop("Numbers of chromosomes between haplotypes are not equal")
    }
    if(identical(boundaries,boundaries2) == F){
      stop("Lengths of chromosomes between haplotypes are not equal")
    }
    digit <- 0
    for(b in 1:length(haplotype1)){
      if(haplotype1[b]==sep){
        digit <- digit + 1
        boundaries[digit] <- boundaries[digit]-(digit-1)
      }
    }
    haplotype1 <- haplotype1[haplotype1 != sep]
    haplotype2 <- haplotype2[haplotype2 != sep]
  } else {
    boundaries.temporary <- c(rep(0,times=length(l.chromosomes)))
    boundaries <- c()
    for(a in 1:(length(l.chromosomes)-1)){
      boundaries[a] <- l.chromosomes[a] + boundaries.temporary[a]
      boundaries.temporary[a+1] <- boundaries[a]
    }
  }
  # number of chiasms on each chromosome
  n.chiasms <- rpois((length(boundaries)+1),lambda = lambda)
  #crossing
  chromosomes1 <- list()
  chromosomes2 <- list()
  if(chiasm.before==T){
    boundaries.plus <- c(0,boundaries,length(haplotype1))
    for(i in 1:(length(boundaries)+1)){
      chromosomes1[[i]] <- haplotype1[(boundaries.plus[i]+1):boundaries.plus[i+1]]
      chromosomes2[[i]] <- haplotype2[(boundaries.plus[i]+1):boundaries.plus[i+1]]
      chiasm.places <- c()
      chiasm.places <- sort(sample(2:length(chromosomes1[[i]]),size = n.chiasms[i]))
      chiasm.places.plus <- c(1,chiasm.places,(length(chromosomes1[[i]])+1))
      for(j in 1:n.chiasms[i]){
        if(j %% 2 != 0){
          chromosomes1.lost <- chromosomes1[[i]][chiasm.places.plus[j]:(chiasm.places.plus[j+1]-1)]
          chromosomes1[[i]][chiasm.places.plus[j]:(chiasm.places.plus[j+1]-1)] <- chromosomes2[[i]][chiasm.places.plus[j]:(chiasm.places.plus[j+1]-1)]
          chromosomes2[[i]][chiasm.places.plus[j]:(chiasm.places.plus[j+1]-1)] <- chromosomes1.lost
        }
      }
    }
  }
  else if(chiasm.before==F){
    boundaries.plus <- c(0,boundaries,length(haplotype1))
    for(i in 1:(length(boundaries)+1)){
      chromosomes1[[i]] <- haplotype1[(boundaries.plus[i]+1):boundaries.plus[i+1]]
      chromosomes2[[i]] <- haplotype2[(boundaries.plus[i]+1):boundaries.plus[i+1]]
      chiasm.places <- c()
      chiasm.places <- sort(sample(1:(length(chromosomes1[[i]])-1),size = n.chiasms[i]))
      chiasm.places.plus <- c(0,chiasm.places,length(chromosomes1[[i]]))
      for(j in 1:n.chiasms[i]){
        if(j %% 2 != 0){
          chromosomes1.lost <- chromosomes1[[i]][(chiasm.places.plus[j]+1):(chiasm.places.plus[j+1])]
          chromosomes1[[i]][(chiasm.places.plus[j]+1):(chiasm.places.plus[j+1])] <- chromosomes2[[i]][(chiasm.places.plus[j]+1):(chiasm.places.plus[j+1])]
          chromosomes2[[i]][(chiasm.places.plus[j]+1):chiasm.places.plus[j+1]] <- chromosomes1.lost
        }
      }
    }
  }
  #tworzenie gamety
  randoms <- sample(1:2,size = (length(boundaries)+1),replace = T)
  gamete <- numeric(0)
  gamete1 <- numeric(0)
  gamete2 <- numeric(0)
  if(sep == ""){
    if(meiotic==F){
      for(k in 1:(length(boundaries)+1)){
        if(randoms[k]==1){
          gamete <- c(gamete,chromosomes1[[k]])
        } else if(randoms[k]==2){
          gamete <- c(gamete,chromosomes2[[k]])
        }
      }
    }
    else if(meiotic==T){
      for(k in 1:(length(boundaries)+1)){
        if(randoms[k]==1){
          gamete1 <- c(gamete1,chromosomes1[[k]])
          gamete2 <- c(gamete2,chromosomes2[[k]])
        } else if(randoms[k]==2){
          gamete1 <- c(gamete1,chromosomes2[[k]])
          gamete2 <- c(gamete2,chromosomes1[[k]])
        }
      }
      gamete <- list(gamete1,gamete2,gamete1,gamete2)
    }
  }
  if(sep != ""){
    if(meiotic==F){
      for(k in 1:(length(boundaries)+1)){
        if(randoms[k]==1){
          gamete <- c(gamete,sep,chromosomes1[[k]])
        } else if(randoms[k]==2){
          gamete <- c(gamete,sep,chromosomes2[[k]])
        }
      }
      gamete <- gamete[2:length(gamete)]
    }
    else if(meiotic==T){
      for(k in 1:(length(boundaries)+1)){
        if(randoms[k]==1){
          gamete1 <- c(gamete1,sep,chromosomes1[[k]])
          gamete2 <- c(gamete2,sep,chromosomes2[[k]])
        } else if(randoms[k]==2){
          gamete1 <- c(gamete1,sep,chromosomes2[[k]])
          gamete2 <- c(gamete2,sep,chromosomes1[[k]])
        }
      }
      gamete1 <- gamete1[2:length(gamete1)]
      gamete2 <- gamete2[2:length(gamete2)]
      gamete <- list(gamete1,gamete2,gamete1,gamete2)
    }
  }
  return(gamete)
}
