options(show.error.locations=TRUE)
# TO CHOOSE!!! -----------------------------------

mating.system <- c(1,2,3,4,5,6)
#mating.system <- c(1)
# "1" - sperm with the lowest number of mutations
# "2" - random mating
# "3" - sperm with the lowest number of mutations in the same loci in comparison to female
# "4" - sperm with the highest fitness (without dominance)
# "5" - sperm with the highest fitness (with dominance included)
# "6" - sperm which will produce offspring with highest fitness
selection.type <- c("constant","variable")
#selection.type <- c("variable")
#constant or variable
dominance.type <- c("constant","half_variable","variable")
#dominance.type <- "constant"
#constant,half_variable (for each s one h),variable (drawn from distribution)
# Functions --------------------------------------------------------------

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
ALLEL.MUT<-function(s2,h2,u){
  a<-s2*(1-2*h2)
  if(a < 0){
    a <- 0
  }
  b<-s2*h2*(1+u)
  c<-(-u)
  if(a == 0){
    q <- (-c)/b
  } else {
    delta<-((b^2)-(4*a*c))
    if(delta==0){
      q<-(-b/(2*a))
    }else if(delta>0){
      q<-(-b-sqrt(delta))/(2*a)
      if(q<0){
        q<-(-b+sqrt(delta))/(2*a) 
      }  
    }
  }
  if(q<0){q <- 0}
  # q wychodzi mniejsze od 1 dla s>0.00013
  q
} #calculation of q - frequency of mutated allel
osobnik.mut<-function(n){
  zmutuja<-samplex(x=loci.sel, size=l.mutacji[n]) 
  for(i in zmutuja){
    if(populacja[n,i]=="mm"){
      populacja[n,i]<-samplex(c("Mm","mM"),size=1)
    }else if (populacja[n,i]=="Mm" | populacja[n,i]=="mM"){
      populacja[n,i]<-"MM"   
    }
  }
  return(populacja)  
} #making mutations
rozklad.s <- function(p = 0.05,loci.number){
  #p - to srednie s
  s <- seq(from=0.01, to=1, by=0.01)
  wynik.s <- (1/p)*exp(-s/p)
  probabilities <- wynik.s/sum(wynik.s)
  samplex(s,size = loci.number, prob = probabilities, replace = T)
}
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
rozklad.h <- function(sztywne = T,srednie.h = 0.36,min.constant = 0.01,max.constant = 100,min.h = 0.01){
  sekwencja <- seq(min.constant,max.constant,by=min.constant) 
  s <- as.numeric(na.omit(selection.coeff))
  g <- c()
  for(k in sekwencja){
    h1 <- (exp(-k*s))/2
    mean <- round(mean(h1),digits = (nchar(min.h)-2))
    if(mean == srednie.h){
      g[k] <- k 
    }
  }
  g <- as.numeric(na.omit(g))
  if(length(g)==0){
    stop("nie znaleziono stalej w rownaniu - zmien wartosci stalej")
  }
  h1 <- numeric(length(g))
  for(k in 1:length(g)){
    h1[k] <- mean((exp(-g[k]*s))/2)
  }
  differences <- abs(srednie.h - h1)
  g.place <- which.min(differences)
  k <- g[g.place]
  if(sztywne == T){
    dominance_coeff <- round((exp(-k*selection.coeff))/2,digits = (nchar(min.h)-2))
  } else {
    dominance_coeff <- c()
    for(i in 1:length(selection.coeff)){
      if(is.na(selection.coeff[i])==T){
        dominance_coeff[i] <- NA
      } else {
        max.h <- round(exp(-k*selection.coeff[i]),digits = (nchar(min.h)-2))
        if(max.h < min.h){
          dominance_coeff[i] <- 0
        } else {
          mozliwosci <- seq(min.h,max.h,by = (10^(-(nchar(min.h)-2))))
          dominance_coeff[i] <- samplex(mozliwosci,size=1)
        }
      }
    }
  }
  dominance_coeff
}

# Given parameters --------------------------------------------------------

N<-10 #size of initial population (even number)
M<-2900*2 #size of haploid genome uder selection (in loci, 2 chromosomes)
L<-2 #average number of chiasms on each chromosome
G<-6000 #overall size of haploid genome (in loci, 200 neutral loci included)
Ul<-0.03 #lethal mutations rate (from Wang's article)
U<-1 #deleterious mutations rate (from Wang's article)
l.pokolen<-50
l.replikatow <- 21
l.plemnikow<-100
l.jaj<-10
s<-0.05 #mean selection coefficient
h<-0.36 #mean dominance coefficient
u<-(U+Ul)/(2*M) #average rate of mutation per locus
no_cores <- 21

# Estimations matrixes formation --------------------------------------------
#zapętlenia
writeLines("","log.txt")
for(mat_sys in mating.system){
  for(s.type in selection.type){
    for(h.type in dominance.type){
      file.name1 <- paste("parametry","mat",mat_sys,"s",s.type,"h",h.type,".txt",sep="")
      parametry<-matrix(nrow=l.pokolen+1, ncol=13) #for each generation across replicates
      colnames(parametry) <- c("mean_viability","max_viability","min_viability","SD_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","alive_offspring","replicates_number")
      row.names(parametry) <- c(0:l.pokolen)
      
      # Ancestral population formation ------------------------------------------
      
      Nmax <- N
      N.fem<-N/2
      N.mal<-N/2
      loci.neutralne<-c(1:200)*30
      wszystkie.loci<-c(1:G)
      loci.sel<-wszystkie.loci[!is.element(wszystkie.loci, loci.neutralne)]
      allele.neutralne<-c("A", "B", "C", "D")
      genotypy.neutralne<-expand.grid(allele.neutralne,allele.neutralne) #creation of all possible combinations
      heterozygoty.neutralne<-subset(genotypy.neutralne, !(Var1==Var2))
      heterozygoty<-paste(heterozygoty.neutralne[,1],heterozygoty.neutralne[,2], sep="")
      
      replicates <- c(1:l.replikatow)
      library(foreach) #introduction to foreach
      library(doParallel)
      #no_cores <- detectCores() - 1
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      parametry.pokolenia <- foreach(rep = replicates) %dopar% {
        if(s.type == "constant"){
          selection.coeff <- rep(s, times = G)
          selection.coeff[loci.neutralne] <- NA
          selection.coeff <- c(selection.coeff,NA)
        } else {
          selection.coeff <- rozklad.s(p=s,loci.number = G)
          selection.coeff[loci.neutralne] <- NA
          selection.coeff <- c(selection.coeff,NA)
        }
        if(h.type == "constant"){
          dominance.coeff <- selection.coeff
          dominance.coeff[!is.na(dominance.coeff)] <- h
        } else if(h.type == "half_variable"){
          dominance.coeff <- rozklad.h(sztywne = T,srednie.h = h)
        } else if(h.type == "variable"){
          dominance.coeff <- rozklad.h(sztywne = F,srednie.h = h)
        }
        populacja<-matrix(nrow=N, ncol=G)
        for(locus in 1:G){
          if(is.element(locus,loci.sel)==T){
            q<-ALLEL.MUT(s2 = selection.coeff[locus],h2 = dominance.coeff[locus],u = u)
            p<-1-q
            cz.MM<-q^2 #frequencies
            cz.mm<-p^2
            cz.Mm<-2*p*q
            populacja[,locus] <- samplex(x=c("MM","mm","Mm", "mM"), size=N, replace=TRUE, prob=c(cz.MM, cz.mm, cz.Mm/2, cz.Mm/2))
          } else {
            string1 <- samplex(allele.neutralne,size = N, replace = T)
            string2 <- samplex(allele.neutralne,size = N, replace = T)
            populacja[,locus] <- paste(string1,string2,sep = "")
          }
        }
        
        # Estimations of ancestral populations ------------------------------------
        
        dostosowanie.populacji <- c() #mean viability
        for(m in 1:nrow(populacja)){
          wierszyk <- populacja[m,1:ncol(populacja)]
          for(n in 1:ncol(populacja)){
            if(wierszyk[n] == "mm"){wierszyk[n] <- 1
            }else if(wierszyk[n] == "MM"){wierszyk[n] <- 1-selection.coeff[n]
            }else if(wierszyk[n] == "Mm" | wierszyk[n]== "mM"){wierszyk[n] <- 1-(dominance.coeff[n]*selection.coeff[n])
            }else {wierszyk[n] <- 1}
          }
          wierszyk <- as.numeric(wierszyk)
          dostosowanie.populacji[m] <- prod(wierszyk)
        }
        mean.viability <- mean(dostosowanie.populacji)
        
        heterozygotycznosc<-c() #mean heterozygosity - number of heterozygotes / number of loci
        l.alleli<-c() #mean number of allels
        fixed.loci<-c(rep(0,times=ncol(populacja))) #number of selective loci with fixed mutation (M)
        segregating.loci<-c(rep(0,times=ncol(populacja))) #number of selective loci with variation (M & m) 
        for(i in loci.sel){
          pierwsza<-substring(populacja[,i],0,1)
          druga<-substring(populacja[,i],2)
          heterozygotycznosc[i]<-sum(pierwsza!=druga)/N
          
          oba<-c(pierwsza,druga)
          l.alleli[i]<-length(unique(oba))
          if(length(unique(oba))==1 & !is.element("m",oba)){fixed.loci[i]<-1}
          if(length(unique(oba))>1 & is.element("M",oba)){segregating.loci[i]<-1}
        }
        mean.heterozygosity <- mean(na.omit(heterozygotycznosc))
        mean.l.alleli<-mean(na.omit(l.alleli))
        mean.fixed_loci<-sum(na.omit(fixed.loci))
        mean.segregating_loci<-sum(na.omit(segregating.loci))
        
        hetero.mut <- c() #mean number of heterozygous mutations per indvidual
        homo.mut <- c() #mean number of homozygous mutations per indvidual
        for(i in 1:nrow(populacja)){
          hetero.mut[i]<-sum(populacja[i,]=="Mm"|populacja[i,]=="mM")
          homo.mut[i]<-sum(populacja[i,]=="MM")
        }
        mean.hetero.mut<-mean(hetero.mut)
        mean.homo.mut<-mean(homo.mut)
        
        offspring_number <- nrow(populacja)
        mutations_per_individual <- (sum(populacja=="Mm"|populacja=="mM")+(sum(populacja=="MM")*2))/offspring_number
        
        wektor.plci <- rep(0,times=Nmax)
        wektor.plci[samplex(c(1:Nmax),size = (Nmax/2))] <- 1
        populacja <- cbind(populacja,wektor.plci)
        populacja <- rbind(populacja,selection.coeff,dominance.coeff,deparse.level = 0)
        
        file.name <- paste("replicate",rep,".txt",sep="")
        write.table(populacja, file=file.name, sep="\t", quote=T)
        wrzutka <- c(mean.viability,mean.heterozygosity,mean.l.alleli,mean.hetero.mut,mean.homo.mut,mean.fixed_loci,mean.segregating_loci,mutations_per_individual,offspring_number)
        wrzutka
      }
      stopCluster(cl)
      
      parametry.pokolenia <- do.call(rbind,parametry.pokolenia)
      colnames(parametry.pokolenia) <- c("mean_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","alive_offspring")
      
      # Mean viability across replicates - for further comparisons --------------
      sr.dostosowanie.populacji <- mean(parametry.pokolenia[,"mean_viability"])
      parametry["0","mean_viability"] <- 1
      parametry["0","max_viability"] <- max((parametry.pokolenia[,"mean_viability"]/sr.dostosowanie.populacji),na.rm = T)
      parametry["0","min_viability"] <- min((parametry.pokolenia[,"mean_viability"]/sr.dostosowanie.populacji),na.rm = T)
      parametry["0","SD_viability"] <- sd((parametry.pokolenia[,"mean_viability"]/sr.dostosowanie.populacji),na.rm = T)
      parametry["0","replicates_number"] <- nrow(parametry.pokolenia)
      for(kol in 2:ncol(parametry.pokolenia)){
        parametry["0",kol+3] <- mean(parametry.pokolenia[,kol])
      }
      write.table(parametry, file=file.name1,sep="\t",quote=T)
      
      # SIMULATION --------------------------------------------------------------
      
      for(pok in 1:l.pokolen){
        if(sum(!is.na(parametry[pok,]))==0){
          break
        }
        
        replicates <- c(1:l.replikatow)
        library(foreach) #introduction to foreach
        library(doParallel)
        #no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores)
        registerDoParallel(cl)
        
        parametry.pokolenia <- foreach(rep = replicates) %dopar% {
          sink("log.txt",append = T)
          cat(paste("pokolenie",pok,"replikat",rep,"combination",mat_sys,s.type,h.type,"\n",sep = " "))
          sink()
          file.name <- paste("replicate",rep,".txt",sep="")
          populacja <- read.table(file.name,sep = "\t",header = T,colClasses = "character")
          selection.coeff <- as.numeric(populacja[(nrow(populacja)-1),])
          dominance.coeff <- as.numeric(populacja[nrow(populacja),])
          populacja <- populacja[1:(nrow(populacja)-2),]
          wrzutka <- c()
          if(nrow(populacja)>=2 & length(unique(populacja[,(G+1)]))==2){ #zabezpieczenie
            N <- nrow(populacja)
            samce <- samplex(which(populacja[,(G+1)]==1))
            samice <- samplex(c(1:N)[!is.element(c(1:N),samce)],replace = F)
            populacja <- populacja[,(1:G)]
            #-------1. Mutations---------
            l.mutacji<-rpois(N, U+Ul)
            for(n in 1:N){
              if(l.mutacji[n]>G){
                l.mutacji[n]<-0
                next
                #populacja<-osobnik.mut(n)
              }else if(l.mutacji[n]==0){
                next
              }else{
                populacja <- osobnik.mut(n)
              }
            }
            #-------2. Mating---------
            chromosomy <- matrix(nrow=2*N, ncol=G) #2 rows per one individual
            for(i in 1:N){
              for(j in 1:G){
                chromosomy[2*i-1,j]<-substring(populacja[i,j],0,1)
                chromosomy[2*i,j]<-substring(populacja[i,j],2)
              }
            }
            plemniki<-matrix(nrow=(length(samce)*l.plemnikow),ncol=G)
            jaja<-matrix(nrow=(length(samice)*l.jaj),ncol=G)
            
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
              plemniki<-cbind(plemniki,ranking.plemnikow)
              if(nrow(plemniki)>=1){
                plemniki <- plemniki[samplex(c(1:nrow(plemniki)), rep=F),]
                plemniki <- plemniki[order(plemniki[,6001],decreasing = F),]
                plemniki <- plemniki[,1:6000]
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
              for(jajo in 1:nrow(jaja)){
                for(plemnik in 1:nrow(plemniki)){
                  odleglosci[jajo,plemnik]<- sum(is.element(jaja.mut[[jajo]],plemniki.mut[[plemnik]]))
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
                nr.plemnika <- which.min(odleglosci[y,])
                odleglosci[,nr.plemnika]<-NA
                gameta2 <- jaja[y,]
                gameta1 <- plemniki[nr.plemnika,]
                for(x in 1:G){
                  populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
                }
              }
            } else if(mat_sys == 4 | mat_sys == 5){
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
                  sperm_fitness <- rep(1,times=G)
                  for(locus in 1:G){
                    if(plemniki[((i-1)*l.plemnikow)+nr,locus]=="M"){
                      if(mat_sys == 4){
                        sperm_fitness[locus] <- 1-selection.coeff[locus]
                      } else if (mat_sys == 5){
                        sperm_fitness[locus] <- 1-(selection.coeff[locus]*dominance.coeff[locus])
                      }
                    }
                  }
                  ranking.plemnikow[((i-1)*l.plemnikow)+nr] <- prod(sperm_fitness)
                }
              }
              plemniki<-cbind(plemniki,ranking.plemnikow)
              if(nrow(plemniki)>=1){
                plemniki <- plemniki[samplex(c(1:nrow(plemniki)), rep=F),]
                plemniki <- plemniki[order(plemniki[,6001],decreasing = T),]
                plemniki <- plemniki[,1:6000]
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
            } else if(mat_sys == 6){
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
              
              #Macierz odleglosci, czyli ile homozygot zmutowanych dalaby kazda z kombinacji jaja z pleminikiem.
              odleglosci<-matrix(nrow=nrow(jaja),ncol=nrow(plemniki))
              effect_homo <- c()
              effect_hetero <- c()
              for(jajo in 1:nrow(jaja)){
                for(plemnik in 1:nrow(plemniki)){
                  place_homo <- jaja.mut[[jajo]][is.element(jaja.mut[[jajo]],plemniki.mut[[plemnik]])]
                  if(length(place_homo)>=1){
                    effect_homo <- 1-selection.coeff[place_homo]
                  }
                  place_hetero <- c(jaja.mut[[jajo]][!is.element(jaja.mut[[jajo]],plemniki.mut[[plemnik]])],plemniki.mut[[plemnik]][!is.element(plemniki.mut[[plemnik]],jaja.mut[[jajo]])])
                  if(length(place_hetero)>=1){
                    effect_hetero <- 1-(selection.coeff[place_hetero]*dominance.coeff[place_hetero])
                  }
                  odleglosci[jajo,plemnik]<- prod(c(effect_homo,effect_hetero))
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
                nr.plemnika <- which.max(odleglosci[y,])
                odleglosci[,nr.plemnika]<-NA
                gameta2 <- jaja[y,]
                gameta1 <- plemniki[nr.plemnika,]
                for(x in 1:G){
                  populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
                }
              }
            }
            
            populacja1 <- cbind(populacja1,c(samplex(c(0,1),size=nrow(populacja1),replace = T)))
            
            # Simulation estimators ---------------------------------------------------
            
            viability <- matrix(nrow=nrow(jaja),ncol = 1)
            for(m in 1:nrow(populacja1)){
              wiersz <- populacja1[m,1:(ncol(populacja1)-1)]
              for(n in 1:(ncol(populacja1)-1)){
                if(wiersz[n] == "mm"){wiersz[n] <- 1
                }else if(wiersz[n] == "MM"){wiersz[n] <- 1-selection.coeff[n]
                }else if(wiersz[n] == "Mm" | wiersz[n]== "mM"){wiersz[n] <- 1-(dominance.coeff[n]*selection.coeff[n])
                }else {wiersz[n] <- 1}
              }
              wiersz <- as.numeric(wiersz)
              viability[m,1] <- prod(wiersz)/sr.dostosowanie.populacji
            }
            survivors <- matrix(nrow=nrow(jaja),ncol=1)
            for (o in 1:nrow(viability)){
              random <- runif(1)
              if(viability[o,1] >= 0.1 & random < viability[o,1]){
                survivors[o,1] <- viability[o,1]
              }
            }
            survivors <- cbind(survivors,populacja1)
            survivors <- survivors[complete.cases(survivors),]
            if(is.null(nrow(survivors))==T){
              survivors <- matrix(survivors,nrow = 1,ncol=(length(survivors)))
            }
            new_population <- survivors
            
            mean.viability<-mean(as.numeric(new_population[,1]))
            
            wektor.plci <- new_population[,(G+2)]
            populacja<-new_population[,2:(G+1)] #for new generation
            if(is.null(nrow(populacja))==T){
              populacja <- matrix(populacja,nrow = 1,ncol=(length(populacja)))
            }
            if(nrow(populacja)>0){
              heterozygotycznosc <- c()
              l.alleli <- c()
              
              fixed.loci<-c(rep(0,times=length(loci.sel)))
              segregating.loci<-c(rep(0,times=length(loci.sel)))
              for(i in loci.sel){
                pierwsza<-substring(populacja[,i],0,1)
                druga<-substring(populacja[,i],2)
                heterozygotycznosc[i]<-sum(pierwsza!=druga)/nrow(populacja)
                
                oba<-c(pierwsza,druga)
                l.alleli[i]<-length(unique(oba))
                if(length(unique(oba))==1 & !is.element("m",oba)){fixed.loci[i]<-1}
                if(length(unique(oba))>1 & is.element("M",oba)){segregating.loci[i]<-1}
              }
              
              mean.heterozygosity <- mean(na.omit(heterozygotycznosc))
              mean.l.alleli<-mean(na.omit(l.alleli))
              mean.fixed_loci<-sum(na.omit(fixed.loci))
              mean.segregating_loci<-sum(na.omit(segregating.loci))
              
              hetero.mut <- c()
              homo.mut <- c()
              if(nrow(populacja)>=1){    
                for(i in 1:nrow(populacja)){
                  hetero.mut[i]<-sum(populacja[i,]=="Mm"|populacja[i,]=="mM")
                  homo.mut[i]<-sum(populacja[i,]=="MM")
                }}
              
              mean.hetero.mut <- mean(hetero.mut)
              mean.homo.mut <- mean(homo.mut)
              offspring_number <- nrow(populacja)
              mutations_per_individual <- (sum(populacja=="Mm"|populacja=="mM")+(sum(populacja=="MM")*2))/offspring_number
              
              males <- which(wektor.plci==1)
              females <- which(wektor.plci==0)
              if(length(males)>=(Nmax/2)){
                males.rows <- samplex(males,size = (Nmax/2))
              } else {males.rows <- males}
              if(length(females)>=(Nmax/2)){
                females.rows <- samplex(females,size = (Nmax/2))
              } else {females.rows <- females}
              populacja <- populacja[c(males.rows,females.rows),]
              if(is.null(nrow(populacja))==T){
                populacja <- matrix(populacja,nrow = 1,ncol=(length(populacja)))
              }
              wektor.plci <- c(rep(1,times=length(males.rows)),rep(0,times=length(females.rows)))
              
              wrzutka <- c(mean.viability,mean.heterozygosity,mean.l.alleli,mean.hetero.mut,mean.homo.mut,mean.fixed_loci,mean.segregating_loci,mutations_per_individual,offspring_number)
              
              populacja <- cbind(populacja,wektor.plci)
              if(nrow(populacja)>1){
                populacja <- populacja[samplex(c(1:nrow(populacja)),rep=F),]
              }
            } else {
              populacja <- matrix(nrow = 0,ncol = (G+1))
            }
              populacja <- rbind(populacja,selection.coeff,dominance.coeff,deparse.level = 0)
              write.table(populacja, file=file.name, sep="\t", quote=T)
              wrzutka
          }
        }
        stopCluster(cl)
        
        parametry.pokolenia<-do.call(rbind,parametry.pokolenia)
        if(length(parametry.pokolenia)>0){
          colnames(parametry.pokolenia) <- c("mean_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","alive_offspring")
          #parametry.pokolenia <- na.omit(parametry.pokolenia)
          parametry[pok+1,"max_viability"] <- max(parametry.pokolenia[,"mean_viability"],na.rm = T)
          parametry[pok+1,"min_viability"] <- min(parametry.pokolenia[,"mean_viability"],na.rm = T)
          parametry[pok+1,"SD_viability"] <- sd(parametry.pokolenia[,"mean_viability"],na.rm = T)
          parametry[pok+1,"replicates_number"] <- nrow(parametry.pokolenia)
          for(kol in 1:ncol(parametry.pokolenia)){
            if(kol == 1){
              parametry[pok+1,kol] <- mean(parametry.pokolenia[,kol])
            } else {
              parametry[pok+1,kol+3] <- mean(parametry.pokolenia[,kol])
            }
          }
        }
        write.table(parametry, file=file.name1,sep="\t",quote=T)
      }
      print("finished")
    }
  }
}

