#BLOSy
samplex <- function(x, size = length(x), replace = FALSE, prob = NULL){
  if(length(x)==1 & replace == F & size > 1){
    stop("Cannot take a sample larger than the population when 'replace = FALSE'")
  } else if(length(x)==1 & replace == F & size == 1){
    outcome <- x
  } else if(length(x)==1 & replace == T){
    outcome <- rep(x,times=size)
  } else {
    outcome <- sample(x, size, replace = replace, prob = prob)
  }
  outcome
}
osobnik.mutacja<-function(n,wektor){
  zmutuja<-samplex(x=loci.sel, size=n) 
  for(i in zmutuja){
    if(wektor[i]=="mm"){
      wektor[i]<-samplex(c("Mm","mM"),size=1)
    }else if (wektor[i]=="Mm" | wektor[i]=="mM"){
      wektor[i]<-"MM"   
    }
  }
  return(wektor)  
} #making mutations
choice.prob <- function(wektor,wspolczynnik,liniowo = T,ideal = F){
  if(length(wektor) < 1){
    stop("wektor musi miec dlugosc co najmniej 1")
  }
  if(wspolczynnik <= 1 & liniowo == F){
    stop("argument _wspolczynnik_ musi byc wiekszy od 1")
  }
  if(wspolczynnik <= 0 & liniowo == T){
    stop("argument _wspolczynnik_ musi byc wiekszy od 0")
  }
  if(ideal == F){
    max.wektor <- max(wektor)
    odwrotnosc <- (max.wektor-wektor)+1
    l.cegielek <- 1
    if(max.wektor>0 & liniowo == T){
      for(i in 1:max.wektor){
        l.cegielek[i+1] <- l.cegielek[i]+wspolczynnik
      }
    }
    if(max.wektor>0 & liniowo == F){
      for(i in 1:max.wektor){
        l.cegielek[i+1] <- l.cegielek[i]*wspolczynnik
      }
    }
    suma.cegielek <- sum(l.cegielek[odwrotnosc])
    cegielka <- 1/suma.cegielek
    prob <- l.cegielek[odwrotnosc]*cegielka
  }
  if(ideal == T){
    min.wektor <- min(wektor)
    prob <- rep(0,times=length(wektor))
    miejsca <- which(wektor == min.wektor)
    prob[miejsca] <- (1/length(miejsca)) 
  }
  prob
}
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
      chiasm.places <- sort(samplex(2:length(chromosomes1[[i]]),size = n.chiasms[i]))
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
      chiasm.places <- sort(samplex(1:(length(chromosomes1[[i]])-1),size = n.chiasms[i]))
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
  randoms <- samplex(1:2,size = (length(boundaries)+1),replace = T)
  gamete <- numeric(0)
  gamete1 <- numeric(0)
  gamete2 <- numeric(0)
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
  return(gamete)
}

mat_sys <- 4
G <- 6000
l.replikatow <- 3
l.jaj <- 10
l.plemnikow <- 100
U<-0.5 #deleterious mutations rate (Simmons, 1977)
loci.neutralne<-c(1:200)*30
wszystkie.loci<-c(1:G)
choice.coeff <- 2
liniowo.value <- T
neutral.allels.nr <- 7
inbreeding.loci.nr <- 5 #how many loci is "mating.system 4" based of
loci.sel<-wszystkie.loci[!is.element(wszystkie.loci, loci.neutralne)]
allele.neutralne<-LETTERS[1:neutral.allels.nr]
mat4.loci <- sort(loci.neutralne[samplex(c(1:length(loci.neutralne)),size = inbreeding.loci.nr,replace = F)])
#I stage--------------------------------------------------------------

for(rep in 1:l.replikatow){
  pop.name1 <- paste("replicate",rep,"samce",sep="")
  pop.name2 <- paste("replicate",rep,"samice",sep="")
  file.name <- paste("replicate",rep,".txt",sep="")
  populacja <- read.table(file.name,sep = "\t",header = T,colClasses = "character")
  selection.coeff <- as.numeric(populacja[(nrow(populacja)-1),])
  dominance.coeff <- as.numeric(populacja[nrow(populacja),])
  populacja <- populacja[1:(nrow(populacja)-2),]
  N <- nrow(populacja)
  samce <- samplex(which(populacja[,(G+1)]==1))
  samice <- samplex(c(1:N)[!is.element(c(1:N),samce)],replace = F)
  samce <- populacja[samce,1:G]
  rownames(samce) <- NULL
  samice <- populacja[samice,1:G]
  rownames(samice) <- NULL
  assign(pop.name1,samce)
  assign(pop.name2,samice)
  
}
rm(populacja)

random.replicates <- samplex(c(1:l.replikatow))
for(rep in random.replicates){
  pop.name2 <- paste("replicate",rep,"samice",sep="")
  samice <- eval(parse(text = pop.name2))
  samica <- samplex(c(1:nrow(samice)),size = 1)
  samica <- as.character(samice[samica,])
  wybor.checking <- F
  double.check <- c()
  while(wybor.checking == F){
    samiec <- c()
    wybor.samca <- c(1:l.replikatow)[-rep]
    wybrany.replikat <- samplex(wybor.samca,size = 1)
    double.check <- sort(unique(c(double.check,wybrany.replikat)))
    pop.name1 <- paste("replicate",wybrany.replikat,"samce",sep="")
    samce <- eval(parse(text = pop.name1))
    if(nrow(samce)>=1){
      samiec.miejsce <- samplex(c(1:nrow(samce)),size = 1)
      samiec <- as.character(samce[samiec.miejsce,])
      samce <- samce[-samiec.miejsce,]
      if(is.null(nrow(samce))==T){
        samce <- matrix(samce,nrow = 1,ncol = length(samce))
      }
      assign(pop.name1,samce)
      wybor.checking <- T
    }
    if(identical(wybor.samca,double.check)==T){
      wybor.checking <- T
    }
  }
  
  if(length(samiec)>0){
    l.mutacji<-rpois(2, U)
    for(sex in 1:2){
      if(l.mutacji[sex]>G){
        l.mutacji[sex]<-0
        next
      }else if(l.mutacji[sex]==0){
        next
      }else if(sex==1){
        samiec <- osobnik.mutacja(n = l.mutacji[sex],samiec)
      }else if(sex==2){
        samica <- osobnik.mutacja(n = l.mutacji[sex],samica)
      }
    }
    #-------2. Mating---------
    chromosomy <- matrix(nrow=4, ncol=G) #2 rows per one individual
    for(sex in 1:2){
      for(j in 1:G){
        if(sex == 1){
          chromosomy[2*sex-1,j]<-substring(samiec[j],0,1)
          chromosomy[2*sex,j]<-substring(samiec[j],2)
        }
        if(sex == 2){
          chromosomy[2*sex-1,j]<-substring(samica[j],0,1)
          chromosomy[2*sex,j]<-substring(samica[j],2)
        }
      }
    }
    plemniki<-matrix(nrow=l.plemnikow,ncol=G)
    jaja<-matrix(nrow=l.jaj,ncol=G)
    samce <- 1
    samice <- 2
    
    if(mat_sys == 1){
      for(i in 1:length(samice)){
        nic2<-samice[i]
        C <- chromosomy[(2*nic2-1),]
        D <- chromosomy[(2*nic2),]
        for(nr in 1:l.jaj){
          jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
          
        }
      }
      if(nrow(jaja)>=1){
        jaja<-jaja[samplex(c(1:nrow(jaja)), rep=F),]
      }
      
      ranking.plemnikow <- c()
      for(i in 1:length(samce)){
        nic<-samce[i]
        A <- chromosomy[(2*nic-1),]
        B <- chromosomy[(2*nic),]
        for(nr in 1:l.plemnikow){
          plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
          ranking.plemnikow[((i-1)*l.plemnikow)+nr] <- sum(plemniki[((i-1)*l.plemnikow)+nr,] == "M")
        }
      }
      #ranking.plemnikow <- choice.prob(ranking.plemnikow,choice.coeff)
      
      if(nrow(jaja)>nrow(plemniki)){
        jaja <- jaja[(1:nrow(plemniki)),]
      }
      
      plemniki<-cbind(plemniki,ranking.plemnikow)
      if(nrow(plemniki)>=1){
        plemniki <- plemniki[samplex(c(1:nrow(plemniki)),rep=F),]
      }
      
      populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
      for(y in 1:nrow(jaja)){
        gameta2 <- jaja[y,]
        ranking.prob <- choice.prob(as.numeric(plemniki[,(G+1)]),choice.coeff,liniowo = liniowo.value)
        miejsce.plemnika <- samplex(c(1:nrow(plemniki)),size=1,prob = ranking.prob)
        gameta1 <- plemniki[miejsce.plemnika,1:G]
        plemniki <- plemniki[-miejsce.plemnika,]
        for(x in 1:G){
          populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
        }
      }
    } else if(mat_sys == 2){
      for(i in 1:length(samice)){
        nic2<-samice[i]
        C <- chromosomy[(2*nic2-1),]
        D <- chromosomy[(2*nic2),]
        for(nr in 1:l.jaj){
          jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
          
        }
      }
      if(nrow(jaja)>=1){
        jaja<-jaja[samplex(c(1:nrow(jaja)), rep=F),]
      }
      
      for(i in 1:length(samce)){
        nic<-samce[i]
        A <- chromosomy[(2*nic-1),]
        B <- chromosomy[(2*nic),]
        for(nr in 1:l.plemnikow){
          plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
        }
      }
      if(nrow(plemniki)>=1){
        plemniki <- plemniki[samplex(c(1:nrow(plemniki)), rep=F),] 
      }
      if(nrow(jaja)>nrow(plemniki)){
        jaja <- jaja[(1:nrow(plemniki)),]
      }
      populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
      for(y in 1:nrow(jaja)){
        gameta2 <- jaja[y,]
        gameta1 <- plemniki[y,]
        for(x in 1:G){
          populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
        }
      }
    } else if(mat_sys == 3){
      plemniki.mut<-list()
      jaja.mut <- list()
      
      for(i in 1:length(samce)){
        nic<-samce[i]
        A <- chromosomy[(2*nic-1),]
        B <- chromosomy[(2*nic),]
        for(nr in 1:l.plemnikow){
          plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
          plemniki.mut[[((i-1)*l.plemnikow)+nr]] <-which(plemniki[((i-1)*l.plemnikow)+nr,]=="M")
        }
      }
      
      for(i in 1:length(samice)){
        nic2<-samice[i]
        C <- chromosomy[(2*nic2-1),]
        D <- chromosomy[(2*nic2),]
        for(nr in 1:l.jaj){
          jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
          jaja.mut[[((i-1)*l.jaj)+nr]] <-which(jaja[((i-1)*l.jaj)+nr,]=="M")
        }
      }
      
      #Macierz odleg?o?ci, czyli ile homozygot zmutowanych da?aby ka?da z kombinacji jaja z pleminikiem.
      odleglosci<-matrix(nrow=nrow(jaja),ncol=nrow(plemniki))
      colnames(odleglosci) <- c(1:ncol(odleglosci))
      for(jajo in 1:nrow(jaja)){
        for(plemnik in 1:nrow(plemniki)){
          odleglosci[jajo,plemnik]<- sum(is.element(jaja.mut[[jajo]],plemniki.mut[[plemnik]]))
          if(jajo == 1){
            colnames(odleglosci)[plemnik] <- paste("plemnik",plemnik,sep = "")
          }
        }
      }
      
      #mieszania
      vec1 <- samplex(c(1:nrow(odleglosci)))
      odleglosci <- odleglosci[vec1,]
      jaja <- jaja[vec1,]
      vec2 <- samplex(c(1:ncol(odleglosci)))
      odleglosci <- odleglosci[,vec2]
      plemniki <- plemniki[vec2,]
      
      if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
        jaja <- jaja[(1:nrow(plemniki)),]
      }
      
      #laczenie gamet w zygoty.
      populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
      for(y in 1:nrow(jaja)){
        ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff,liniowo = liniowo.value)
        nr.kol <- samplex(c(1:ncol(odleglosci)),size=1,prob = ranking.prob)
        name.kol <- colnames(odleglosci)[nr.kol]
        nr.plemnika <- as.numeric(substring(name.kol,8))
        odleglosci <- odleglosci[,colnames(odleglosci) != name.kol]
        gameta2 <- jaja[y,]
        gameta1 <- plemniki[nr.plemnika,]
        for(x in 1:G){
          populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
        }
      }
    } else if(mat_sys == 4){
      plemniki.neutral<-list()
      jaja.neutral <- list()
      
      for(i in 1:length(samce)){
        nic<-samce[i]
        A <- chromosomy[(2*nic-1),]
        B <- chromosomy[(2*nic),]
        for(nr in 1:l.plemnikow){
          plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
          plemniki.neutral[[((i-1)*l.plemnikow)+nr]] <- plemniki[((i-1)*l.plemnikow)+nr,mat4.loci]
        }
      }
      
      for(i in 1:length(samice)){
        nic2<-samice[i]
        C <- chromosomy[(2*nic2-1),]
        D <- chromosomy[(2*nic2),]
        for(nr in 1:l.jaj){
          jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
          jaja.neutral[[((i-1)*l.jaj)+nr]] <- jaja[((i-1)*l.jaj)+nr,mat4.loci]
        }
      }
      
      #Macierz odleg?o?ci, czyli ile homozygot zmutowanych da?aby ka?da z kombinacji jaja z pleminikiem.
      odleglosci<-matrix(nrow=nrow(jaja),ncol=nrow(plemniki))
      colnames(odleglosci) <- c(1:ncol(odleglosci))
      for(jajo in 1:nrow(jaja)){
        for(plemnik in 1:nrow(plemniki)){
          odleglosci[jajo,plemnik]<- sum(jaja.neutral[[jajo]] == plemniki.neutral[[plemnik]])
          if(jajo == 1){
            colnames(odleglosci)[plemnik] <- paste("plemnik",plemnik,sep = "")
          }
        }
      }
      
      #mieszania
      vec1 <- samplex(c(1:nrow(odleglosci)))
      odleglosci <- odleglosci[vec1,]
      jaja <- jaja[vec1,]
      vec2 <- samplex(c(1:ncol(odleglosci)))
      odleglosci <- odleglosci[,vec2]
      plemniki <- plemniki[vec2,]
      
      if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
        jaja <- jaja[(1:nrow(plemniki)),]
      }
      
      #laczenie gamet w zygoty.
      populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
      for(y in 1:nrow(jaja)){
        ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff,liniowo = liniowo.value)
        nr.kol <- samplex(c(1:ncol(odleglosci)),size=1,prob = ranking.prob)
        name.kol <- colnames(odleglosci)[nr.kol]
        nr.plemnika <- as.numeric(substring(name.kol,8))
        odleglosci <- odleglosci[,colnames(odleglosci) != name.kol]
        gameta2 <- jaja[y,]
        gameta1 <- plemniki[nr.plemnika,]
        for(x in 1:G){
          populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
        }
      }
    }
    populacja1 <- cbind(populacja1,c(samplex(c(0,1),size=nrow(populacja1),replace = T)))
    file.name4 <- paste("BLOS",rep,".txt",sep = "")
    write.table(populacja1, file=file.name4, sep="\t", quote=T)
  }
}
rm(list=ls(pattern = "^replicate"))
rm(chromosomy,jaja,plemniki,populacja1,odleglosci)

#II stage--------------------------------------------------------------------

for(rep in 1:l.replikatow){
  pop.name1 <- paste("BLOS",rep,"samce",sep="")
  pop.name2 <- paste("BLOS",rep,"samice",sep="")
  file.name <- paste("BLOS",rep,".txt",sep="")
  populacja <- read.table(file.name,sep = "\t",header = T,colClasses = "character")
  N <- nrow(populacja)
  samce <- samplex(which(populacja[,(G+1)]==1))
  samice <- samplex(c(1:N)[!is.element(c(1:N),samce)],replace = F)
  samce <- populacja[samce,1:G]
  rownames(samce) <- NULL
  samice <- populacja[samice,1:G]
  rownames(samice) <- NULL
  assign(pop.name1,samce)
  assign(pop.name2,samice)
}
rm(populacja)

random.replicates <- samplex(c(1:l.replikatow))
breeding <- c("out","in")
for(rep in random.replicates){
  pop.name2 <- paste("BLOS",rep,"samice",sep="")
  samice <- eval(parse(text = pop.name2))
  if(nrow(samice)>=2){
    for(bred in breeding){
      pop.name2 <- paste("BLOS",rep,"samice",sep="")
      samice <- eval(parse(text = pop.name2))
      samica.miejsce <- samplex(c(1:nrow(samice)),size = 1)
      samica <- as.character(samice[samica.miejsce,])
      samice <- samice[-samica.miejsce,]
      if(is.null(nrow(samice))==T){
        samice <- matrix(samice,nrow = 1,ncol = length(samice))
      }
      assign(pop.name2,samice)
      
      if(bred == "out"){
        wybor.checking <- F
        double.check <- c()
        while(wybor.checking == F){
          samiec <- c()
          wybor.samca <- c(1:l.replikatow)[-rep]
          wybrany.replikat <- samplex(wybor.samca,size = 1)
          double.check <- sort(unique(c(double.check,wybrany.replikat)))
          pop.name1 <- paste("BLOS",wybrany.replikat,"samce",sep="")
          samce <- eval(parse(text = pop.name1))
          if(nrow(samce)>=1){
            samiec.miejsce <- samplex(c(1:nrow(samce)),size = 1)
            samiec <- as.character(samce[samiec.miejsce,])
            samce <- samce[-samiec.miejsce,]
            if(is.null(nrow(samce))==T){
              samce <- matrix(samce,nrow = 1,ncol = length(samce))
            }
            assign(pop.name1,samce)
            wybor.checking <- T
          }
          if(identical(wybor.samca,double.check)==T){
            wybor.checking <- T
          }
        }
      }
      
      if(bred == "in"){
        samiec <- c()
        pop.name1 <- paste("BLOS",rep,"samce",sep="")
        samce <- eval(parse(text = pop.name1))
        if(nrow(samce)>=1){
          samiec.miejsce <- samplex(c(1:nrow(samce)),size = 1)
          samiec <- as.character(samce[samiec.miejsce,])
          samce <- samce[-samiec.miejsce,]
          if(is.null(nrow(samce))==T){
            samce <- matrix(samce,nrow = 1,ncol = length(samce))
          }
          assign(pop.name1,samce)
        }
      }

      if(length(samiec)>0){
        l.mutacji<-rpois(2, U)
        for(sex in 1:2){
          if(l.mutacji[sex]>G){
            l.mutacji[sex]<-0
            next
          }else if(l.mutacji[sex]==0){
            next
          }else if(sex==1){
            samiec <- osobnik.mutacja(n = l.mutacji[sex],samiec)
          }else if(sex==2){
            samica <- osobnik.mutacja(n = l.mutacji[sex],samica)
          }
        }
        #-------2. Mating---------
        chromosomy <- matrix(nrow=4, ncol=G) #2 rows per one individual
        for(sex in 1:2){
          for(j in 1:G){
            if(sex == 1){
              chromosomy[2*sex-1,j]<-substring(samiec[j],0,1)
              chromosomy[2*sex,j]<-substring(samiec[j],2)
            }
            if(sex == 2){
              chromosomy[2*sex-1,j]<-substring(samica[j],0,1)
              chromosomy[2*sex,j]<-substring(samica[j],2)
            }
          }
        }
        plemniki<-matrix(nrow=l.plemnikow,ncol=G)
        jaja<-matrix(nrow=l.jaj,ncol=G)
        samce <- 1
        samice <- 2
        
        if(mat_sys == 1){
          for(i in 1:length(samice)){
            nic2<-samice[i]
            C <- chromosomy[(2*nic2-1),]
            D <- chromosomy[(2*nic2),]
            for(nr in 1:l.jaj){
              jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
              
            }
          }
          if(nrow(jaja)>=1){
            jaja<-jaja[samplex(c(1:nrow(jaja)), rep=F),]
          }
          
          ranking.plemnikow <- c()
          for(i in 1:length(samce)){
            nic<-samce[i]
            A <- chromosomy[(2*nic-1),]
            B <- chromosomy[(2*nic),]
            for(nr in 1:l.plemnikow){
              plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
              ranking.plemnikow[((i-1)*l.plemnikow)+nr] <- sum(plemniki[((i-1)*l.plemnikow)+nr,] == "M")
            }
          }
          #ranking.plemnikow <- choice.prob(ranking.plemnikow,choice.coeff)
          
          if(nrow(jaja)>nrow(plemniki)){
            jaja <- jaja[(1:nrow(plemniki)),]
          }
          
          plemniki<-cbind(plemniki,ranking.plemnikow)
          if(nrow(plemniki)>=1){
            plemniki <- plemniki[samplex(c(1:nrow(plemniki)),rep=F),]
          }
          
          populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
          for(y in 1:nrow(jaja)){
            gameta2 <- jaja[y,]
            ranking.prob <- choice.prob(as.numeric(plemniki[,(G+1)]),choice.coeff,liniowo = liniowo.value)
            miejsce.plemnika <- samplex(c(1:nrow(plemniki)),size=1,prob = ranking.prob)
            gameta1 <- plemniki[miejsce.plemnika,1:G]
            plemniki <- plemniki[-miejsce.plemnika,]
            for(x in 1:G){
              populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
            }
          }
        } else if(mat_sys == 2){
          for(i in 1:length(samice)){
            nic2<-samice[i]
            C <- chromosomy[(2*nic2-1),]
            D <- chromosomy[(2*nic2),]
            for(nr in 1:l.jaj){
              jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
              
            }
          }
          if(nrow(jaja)>=1){
            jaja<-jaja[samplex(c(1:nrow(jaja)), rep=F),]
          }
          
          for(i in 1:length(samce)){
            nic<-samce[i]
            A <- chromosomy[(2*nic-1),]
            B <- chromosomy[(2*nic),]
            for(nr in 1:l.plemnikow){
              plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
            }
          }
          if(nrow(plemniki)>=1){
            plemniki <- plemniki[samplex(c(1:nrow(plemniki)), rep=F),] 
          }
          if(nrow(jaja)>nrow(plemniki)){
            jaja <- jaja[(1:nrow(plemniki)),]
          }
          populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
          for(y in 1:nrow(jaja)){
            gameta2 <- jaja[y,]
            gameta1 <- plemniki[y,]
            for(x in 1:G){
              populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
            }
          }
        } else if(mat_sys == 3){
          plemniki.mut<-list()
          jaja.mut <- list()
          
          for(i in 1:length(samce)){
            nic<-samce[i]
            A <- chromosomy[(2*nic-1),]
            B <- chromosomy[(2*nic),]
            for(nr in 1:l.plemnikow){
              plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
              plemniki.mut[[((i-1)*l.plemnikow)+nr]] <-which(plemniki[((i-1)*l.plemnikow)+nr,]=="M")
            }
          }
          
          for(i in 1:length(samice)){
            nic2<-samice[i]
            C <- chromosomy[(2*nic2-1),]
            D <- chromosomy[(2*nic2),]
            for(nr in 1:l.jaj){
              jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
              jaja.mut[[((i-1)*l.jaj)+nr]] <-which(jaja[((i-1)*l.jaj)+nr,]=="M")
            }
          }
          
          #Macierz odleg?o?ci, czyli ile homozygot zmutowanych da?aby ka?da z kombinacji jaja z pleminikiem.
          odleglosci<-matrix(nrow=nrow(jaja),ncol=nrow(plemniki))
          colnames(odleglosci) <- c(1:ncol(odleglosci))
          for(jajo in 1:nrow(jaja)){
            for(plemnik in 1:nrow(plemniki)){
              odleglosci[jajo,plemnik]<- sum(is.element(jaja.mut[[jajo]],plemniki.mut[[plemnik]]))
              if(jajo == 1){
                colnames(odleglosci)[plemnik] <- paste("plemnik",plemnik,sep = "")
              }
            }
          }
          
          #mieszania
          vec1 <- samplex(c(1:nrow(odleglosci)))
          odleglosci <- odleglosci[vec1,]
          jaja <- jaja[vec1,]
          vec2 <- samplex(c(1:ncol(odleglosci)))
          odleglosci <- odleglosci[,vec2]
          plemniki <- plemniki[vec2,]
          
          if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
            jaja <- jaja[(1:nrow(plemniki)),]
          }
          
          #laczenie gamet w zygoty.
          populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
          for(y in 1:nrow(jaja)){
            ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff,liniowo = liniowo.value)
            nr.kol <- samplex(c(1:ncol(odleglosci)),size=1,prob = ranking.prob)
            name.kol <- colnames(odleglosci)[nr.kol]
            nr.plemnika <- as.numeric(substring(name.kol,8))
            odleglosci <- odleglosci[,colnames(odleglosci) != name.kol]
            gameta2 <- jaja[y,]
            gameta1 <- plemniki[nr.plemnika,]
            for(x in 1:G){
              populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
            }
          }
        } else if(mat_sys == 4){
          plemniki.neutral<-list()
          jaja.neutral <- list()
          
          for(i in 1:length(samce)){
            nic<-samce[i]
            A <- chromosomy[(2*nic-1),]
            B <- chromosomy[(2*nic),]
            for(nr in 1:l.plemnikow){
              plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = c((G/2),(G/2)))
              plemniki.neutral[[((i-1)*l.plemnikow)+nr]] <- plemniki[((i-1)*l.plemnikow)+nr,mat4.loci]
            }
          }
          
          for(i in 1:length(samice)){
            nic2<-samice[i]
            C <- chromosomy[(2*nic2-1),]
            D <- chromosomy[(2*nic2),]
            for(nr in 1:l.jaj){
              jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = c((G/2),(G/2)))
              jaja.neutral[[((i-1)*l.jaj)+nr]] <- jaja[((i-1)*l.jaj)+nr,mat4.loci]
            }
          }
          
          #Macierz odleg?o?ci, czyli ile homozygot zmutowanych da?aby ka?da z kombinacji jaja z pleminikiem.
          odleglosci<-matrix(nrow=nrow(jaja),ncol=nrow(plemniki))
          colnames(odleglosci) <- c(1:ncol(odleglosci))
          for(jajo in 1:nrow(jaja)){
            for(plemnik in 1:nrow(plemniki)){
              odleglosci[jajo,plemnik]<- sum(jaja.neutral[[jajo]] == plemniki.neutral[[plemnik]])
              if(jajo == 1){
                colnames(odleglosci)[plemnik] <- paste("plemnik",plemnik,sep = "")
              }
            }
          }
          
          #mieszania
          vec1 <- samplex(c(1:nrow(odleglosci)))
          odleglosci <- odleglosci[vec1,]
          jaja <- jaja[vec1,]
          vec2 <- samplex(c(1:ncol(odleglosci)))
          odleglosci <- odleglosci[,vec2]
          plemniki <- plemniki[vec2,]
          
          if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
            jaja <- jaja[(1:nrow(plemniki)),]
          }
          
          #laczenie gamet w zygoty.
          populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
          for(y in 1:nrow(jaja)){
            ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff,liniowo = liniowo.value)
            nr.kol <- samplex(c(1:ncol(odleglosci)),size=1,prob = ranking.prob)
            name.kol <- colnames(odleglosci)[nr.kol]
            nr.plemnika <- as.numeric(substring(name.kol,8))
            odleglosci <- odleglosci[,colnames(odleglosci) != name.kol]
            gameta2 <- jaja[y,]
            gameta1 <- plemniki[nr.plemnika,]
            for(x in 1:G){
              populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
            }
          }
        }
        populacja1 <- cbind(populacja1,c(samplex(c(0,1),size=nrow(populacja1),replace = T)))
        file.name5 <- paste("BLOS",bred,rep,".txt",sep = "")
        write.table(populacja1, file=file.name5, sep="\t", quote=T)
      }
    }
  }
}
rm(list=ls(pattern = "^BLOS"))
rm(chromosomy,jaja,plemniki,populacja1,odleglosci)





  
