#Symulacja z tworzeniem inicjalnej populacji
#Usuwanie mutacji z małych populacji
#Aktualna na dzień 18.03.2018

options(show.error.locations=TRUE)
# TO CHOOSE!!! -----------------------------------

#mating.system <- c(1,2,3,4)
mating.system <- c(3)
# "1" - sperm with the lowest number of mutations
# "2" - random mating
# "3" - sperm with the lowest number of mutations in the same loci in comparison to female
# "4" - sperm with different allels in set of neutral allels
#selection.type <- c("constant","variable")
selection.type <- c("constant")
#constant or variable
#dominance.type <- c("constant","half_variable","variable")
dominance.type <- "constant"
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
      if(q<0 | q>1){
        q<-(-b+sqrt(delta))/(2*a) 
      }  
    }
  }
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
  s <- seq(from=0.001, to=1, by=0.001)
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
choice.prob <- function(wektor,wspolczynnik,ideal = F){
  if(length(wektor) < 1){
    stop("wektor musi miec dlugosc co najmniej 1")
  }
  if(wspolczynnik <= 1){
    stop("argument _wspolczynnik_ musi byc wiekszy od 1")
  }
  if(ideal == F){
    max.wektor <- max(wektor)
    odwrotnosc <- (max.wektor-wektor)+1
    l.cegielek <- 1
    if(max.wektor>0){
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

# Given parameters --------------------------------------------------------

N<-10 #size of populations (even number)
M<-2900*2 #size of haploid genome uder selection (in loci, 2 chromosomes)
L<-2 #average number of chiasms on each chromosome
G<-6000 #overall size of haploid genome (in loci, 200 neutral loci included)
#Ul<-0.03 #lethal mutations rate (from Wang's article)
U<-0.5 #deleterious mutations rate (Simmons, 1977)
l.pokolen <- 50
l.replikatow <- 21
l.plemnikow<-100
l.jaj<-10
s<-0.05 #mean selection coefficient (Halligan, Keightley, 2009)
h<-0.36 #mean dominance coefficient (Garcia-Dorado,Caballero, 2001)
u<-U/(2*M) #average rate of mutation per locus
no_cores <- l.replikatow
initial.pop <- 10000
#choice.coeff <- 8 #how many times probability of choice decreases with 1 additional mutation (must be larger than 1)
choice.coeff.wektor <- c(1.1,1.15,1.2,1.3,1.35,1.4)
neutral.allels.nr <- 7
inbreeding.loci.nr <- 5 #how many loci is "mating.system 4" based of

# Estimations matrixes formation --------------------------------------------
#zapętlenia
writeLines("","log.txt")
for(mat_sys in mating.system){
  for(s.type in selection.type){
    for(h.type in dominance.type){
      for(choice.coeff in choice.coeff.wektor){
        if(!(s.type == "constant" & h.type == "half_variable")){
          file.name1 <- paste("parametry","mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
          parametry<-matrix(nrow=l.pokolen+1, ncol=14) #for each generation across replicates
          colnames(parametry) <- c("mean_viability","max_viability","min_viability","SD_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","lethal_eq","alive_offspring","replicates_number")
          row.names(parametry) <- c(0:l.pokolen)
          
          # Ancestral population formation ------------------------------------------
          
          Nmax <- N
          N.fem<-N/2
          N.mal<-N/2
          loci.neutralne<-c(1:200)*30
          wszystkie.loci<-c(1:G)
          loci.sel<-wszystkie.loci[!is.element(wszystkie.loci, loci.neutralne)]
          allele.neutralne<-LETTERS[1:neutral.allels.nr]
          mat4.loci <- sort(loci.neutralne[samplex(c(1:length(loci.neutralne)),size = inbreeding.loci.nr,replace = F)])
          
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
          
          populacja<-matrix(nrow=initial.pop, ncol=G+1)
          for(locus in 1:G){
            if(is.element(locus,loci.sel)==T){
              q<-ALLEL.MUT(s2 = selection.coeff[locus],h2 = dominance.coeff[locus],u = u)
              p<-1-q
              cz.MM<-q^2 #frequencies
              cz.mm<-p^2
              cz.Mm<-2*p*q
              populacja[,locus] <- samplex(x=c("MM","mm","Mm", "mM"), size=initial.pop, replace=TRUE, prob=c(cz.MM, cz.mm, cz.Mm/2, cz.Mm/2))
            } else {
              string1 <- samplex(allele.neutralne,size = initial.pop, replace = T)
              string2 <- samplex(allele.neutralne,size = initial.pop, replace = T)
              populacja[,locus] <- paste(string1,string2,sep = "")
            }
          }
          
          library(foreach) #introduction to foreach
          library(doParallel)
          cl <- makeCluster(no_cores)
          registerDoParallel(cl)
          init.pop <- c(1:initial.pop)
          sr.dostosowanie.populacji <- foreach(init = init.pop) %dopar% {
            homo.effect <- c()
            hetero.effect <- c()
            homo.place <- which(populacja[init,]=="MM")
            if(length(homo.place)>=1){
              homo.effect <- 1-selection.coeff[homo.place]
            }
            hetero.place <- which(populacja[init,]=="mM" | populacja[init,]=="Mm")
            if(length(hetero.place)>=1){
              hetero.effect <- 1-(selection.coeff[hetero.place]*dominance.coeff[hetero.place])
            }
            wynik <- prod(c(homo.effect,hetero.effect))
            wynik
          }
          stopCluster(cl)
          sr.dostosowanie.populacji <- as.vector(do.call(cbind,sr.dostosowanie.populacji))
          sr.dostosowanie.populacji <- mean(sr.dostosowanie.populacji)
          
          populacja[,G+1] <- samplex(c(0,1),size = initial.pop,replace = T)
          for(i in 1:l.replikatow){
            males <- which(populacja[,G+1]==1)
            females <- which(populacja[,G+1]==0)
            if(length(males)>=(Nmax/2)){
              males.rows <- samplex(males,size = (Nmax/2))
            } else {males.rows <- males}
            if(length(females)>=(Nmax/2)){
              females.rows <- samplex(females,size = (Nmax/2))
            } else {females.rows <- females}
            rows <- c(males.rows,females.rows)
            pop <- populacja[rows,]
            pop <- rbind(pop,selection.coeff,dominance.coeff,deparse.level = 0)
            file.name <- paste("replicate",i,"mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
            write.table(pop, file=file.name, sep="\t", quote=T)
            populacja <- populacja[-rows,]
          }
          rm(populacja)
          # Estimations of ancestral populations ------------------------------------
          
          replicates <- c(1:l.replikatow)
          library(foreach) #introduction to foreach
          library(doParallel)
          #no_cores <- detectCores() - 1
          cl <- makeCluster(no_cores)
          registerDoParallel(cl)
          
          parametry.pokolenia <- foreach(rep = replicates) %dopar% {
            file.name <- paste("replicate",rep,"mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
            file.name3 <- paste("plot","rep",rep,"mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
            populacja <- read.table(file.name,sep = "\t",header = T,colClasses = "character")
            selection.coeff <- as.numeric(populacja[(nrow(populacja)-1),])
            dominance.coeff <- as.numeric(populacja[nrow(populacja),])
            wektor.plci <- populacja[,ncol(populacja)]
            samce <- which(wektor.plci==1)
            samice <- which(wektor.plci==0)
            populacja <- populacja[1:(nrow(populacja)-2),1:(ncol(populacja)-1)]
            
            x.plot.homo.1 <- 0
            y.plot.homo.1 <- 0
            x.plot.hetero.1 <- 0
            y.plot.hetero.1 <- 0
            x.plot.homo.0 <- 0
            y.plot.homo.0 <- 0
            x.plot.hetero.0 <- 0
            y.plot.hetero.0 <- 0
            invisible(file.remove(file.name3))
            dostosowanie.populacji <- c() #mean viability
            for(m in 1:nrow(populacja)){
              wierszyk <- populacja[m,1:ncol(populacja)]
              for(n in 1:ncol(populacja)){
                if(wierszyk[n] == "mm"){wierszyk[n] <- 1
                }else if(wierszyk[n] == "MM"){
                  wierszyk[n] <- 1-selection.coeff[n]
                  if(is.element(m,samce)==T){
                    x.plot.homo.1[length(x.plot.homo.1)+1] <- n
                    y.plot.homo.1[length(y.plot.homo.1)+1] <- m
                  } else if(is.element(m,samice)==T){
                    x.plot.homo.0[length(x.plot.homo.0)+1] <- n
                    y.plot.homo.0[length(y.plot.homo.0)+1] <- m
                  }
                }else if(wierszyk[n] == "Mm" | wierszyk[n]== "mM"){
                  wierszyk[n] <- 1-(dominance.coeff[n]*selection.coeff[n])
                  if(is.element(m,samce)==T){
                    x.plot.hetero.1[length(x.plot.hetero.1)+1] <- n
                    y.plot.hetero.1[length(y.plot.hetero.1)+1] <- m 
                  } else if(is.element(m,samice)==T){
                    x.plot.hetero.0[length(x.plot.hetero.0)+1] <- n
                    y.plot.hetero.0[length(y.plot.hetero.0)+1] <- m 
                  }
                }else {wierszyk[n] <- 1}
              }
              wierszyk <- as.numeric(wierszyk)
              dostosowanie.populacji[m] <- prod(wierszyk)
            }
            mean.viability <- mean(dostosowanie.populacji)/sr.dostosowanie.populacji
            coordinates <- list(x.plot.homo.1,y.plot.homo.1,x.plot.hetero.1,y.plot.hetero.1,x.plot.homo.0,y.plot.homo.0,x.plot.hetero.0,y.plot.hetero.0)
            coordinates.length <- c()
            for(cor in 1:length(coordinates)){
              coordinates.length[cor] <- length(coordinates[[cor]]) 
            }
            coordinates.length <- max(coordinates.length)
            for(cor in 1:length(coordinates)){
              roznica <- coordinates.length-length(coordinates[[cor]])
              coordinates[[cor]] <- c(coordinates[[cor]],rep(NA,times=roznica)) 
            }
            coordinates <- do.call(rbind,coordinates)
            write.table(coordinates, file=file.name3,sep="\t",quote=T)
            
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
            lethal <- c()
            for(i in 1:nrow(populacja)){
              hetero.mut[i]<-sum(populacja[i,]=="Mm"|populacja[i,]=="mM")
              homo.mut[i]<-sum(populacja[i,]=="MM")
              lethal.hetero <- which(populacja[i,]=="Mm" | populacja[i,]=="mM")
              lethal.homo <- which(populacja[i,]=="MM")
              lethal[i] <- sum(selection.coeff[lethal.hetero])+(sum(selection.coeff[lethal.homo])*2)
            }
            mean.hetero.mut<-mean(hetero.mut)
            mean.homo.mut<-mean(homo.mut)
            lethal_eq <- mean(lethal)
            
            offspring_number <- nrow(populacja)
            mutations_per_individual <- (sum(populacja=="Mm"|populacja=="mM")+(sum(populacja=="MM")*2))/offspring_number
            
            wrzutka <- c(mean.viability,mean.heterozygosity,mean.l.alleli,mean.hetero.mut,mean.homo.mut,mean.fixed_loci,mean.segregating_loci,mutations_per_individual,lethal_eq,offspring_number)
            wrzutka
          }
          stopCluster(cl)
          
          parametry.pokolenia <- do.call(rbind,parametry.pokolenia)
          colnames(parametry.pokolenia) <- c("mean_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","lethal_eq","alive_offspring")
          
          # Mean viability across replicates - for further comparisons --------------
          #sr.dostosowanie.populacji <- mean(parametry.pokolenia[,"mean_viability"])
          parametry["0","mean_viability"] <- mean(parametry.pokolenia[,"mean_viability"])
          parametry["0","max_viability"] <- max(parametry.pokolenia[,"mean_viability"],na.rm = T)
          parametry["0","min_viability"] <- min(parametry.pokolenia[,"mean_viability"],na.rm = T)
          parametry["0","SD_viability"] <- sd(parametry.pokolenia[,"mean_viability"],na.rm = T)
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
              cat(paste("pokolenie",pok,"replikat",rep,"combination",mat_sys,s.type,h.type,choice.coeff,"\n",sep = " "))
              sink()
              file.name <- paste("replicate",rep,"mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
              file.name3 <- paste("plot","rep",rep,"mat",mat_sys,"s",s.type,"h",h.type,"_choice",choice.coeff,".txt",sep="")
              file.name4 <- paste("choice_range",choice.coeff,"rep",rep,"mat",mat_sys,"s",s.type,"h",h.type,".txt",sep="")
              if(pok == 1){
                writeLines("",file.name4)
              }
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
                l.mutacji<-rpois(N, U)
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
                    ranking.prob <- choice.prob(as.numeric(plemniki[,(G+1)]),choice.coeff)
                    if(y == 1){
                      Wincenty_Kadlubek <- as.numeric(plemniki[,(G+1)])
                    }
                    miejsce.plemnika <- samplex(c(1:nrow(plemniki)),size=1,prob = ranking.prob)
                    gameta1 <- plemniki[miejsce.plemnika,1:G]
                    plemniki <- plemniki[-miejsce.plemnika,]
                    for(x in 1:G){
                      populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
                    }
                  }
                  write(Wincenty_Kadlubek,file = file.name4,ncolumns = 500,append = T)
                  
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
                  
                  if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
                    jaja <- jaja[(1:nrow(plemniki)),]
                  }
                  
                  #laczenie gamet w zygoty.
                  populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
                  for(y in 1:nrow(jaja)){
                    if(y == 1){
                      Wincenty_Kadlubek <- odleglosci[y,]
                    }
                    ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff)
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
                  write(Wincenty_Kadlubek,file = file.name4,ncolumns = 500,append = T)
                  
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
                  
                  if(nrow(jaja)>nrow(plemniki)){ #zabezpieczenie
                    jaja <- jaja[(1:nrow(plemniki)),]
                  }
                  
                  #laczenie gamet w zygoty.
                  populacja1 <- matrix(nrow=nrow(jaja),ncol=G)
                  for(y in 1:nrow(jaja)){
                    ranking.prob <- choice.prob(as.numeric(odleglosci[y,]),choice.coeff)
                    if(y == 1){
                      Wincenty_Kadlubek <- odleglosci[y,]
                    }
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
                  write(Wincenty_Kadlubek,file = file.name4,ncolumns = 500,append = T)
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
                  if(random < viability[o,1]){
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
                  lethal <- c()
                  if(nrow(populacja)>=1){    
                    for(i in 1:nrow(populacja)){
                      hetero.mut[i]<-sum(populacja[i,]=="Mm"|populacja[i,]=="mM")
                      homo.mut[i]<-sum(populacja[i,]=="MM")
                      lethal.hetero <- which(populacja[i,]=="Mm" | populacja[i,]=="mM")
                      lethal.homo <- which(populacja[i,]=="MM")
                      lethal[i] <- sum(selection.coeff[lethal.hetero])+(sum(selection.coeff[lethal.homo])*2)
                    }}
                  
                  mean.hetero.mut <- mean(hetero.mut)
                  mean.homo.mut <- mean(homo.mut)
                  lethal_eq <- mean(lethal)
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
                  
                  wrzutka <- c(mean.viability,mean.heterozygosity,mean.l.alleli,mean.hetero.mut,mean.homo.mut,mean.fixed_loci,mean.segregating_loci,mutations_per_individual,lethal_eq,offspring_number)
                  
                  coordinatesTable <- read.table(file.name3,sep = "\t",header = T,colClasses = "character")
                  invisible(file.remove(file.name3))
                  x.plot.homo.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[1,])))
                  y.plot.homo.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[2,])))
                  x.plot.hetero.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[3,])))
                  y.plot.hetero.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[4,])))
                  x.plot.homo.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[5,])))
                  y.plot.homo.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[6,])))
                  x.plot.hetero.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[7,])))
                  y.plot.hetero.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[8,])))
                  samce <- c(1:length(males.rows))
                  samice <- c((length(males.rows)+1):length(wektor.plci))
                  for(m in 1:nrow(populacja)){
                    wierszyk <- populacja[m,1:ncol(populacja)]
                    for(n in 1:ncol(populacja)){
                      if(wierszyk[n] == "MM"){
                        if(is.element(m,samce)==T){
                          x.plot.homo.1[length(x.plot.homo.1)+1] <- n
                          y.plot.homo.1[length(y.plot.homo.1)+1] <- (pok*Nmax)+m
                        } else if(is.element(m,samice)==T){
                          x.plot.homo.0[length(x.plot.homo.0)+1] <- n
                          y.plot.homo.0[length(y.plot.homo.0)+1] <- (pok*Nmax)+m
                        }
                      }else if(wierszyk[n] == "Mm" | wierszyk[n]== "mM"){
                        if(is.element(m,samce)==T){
                          x.plot.hetero.1[length(x.plot.hetero.1)+1] <- n
                          y.plot.hetero.1[length(y.plot.hetero.1)+1] <- (pok*Nmax)+m 
                        } else if(is.element(m,samice)==T){
                          x.plot.hetero.0[length(x.plot.hetero.0)+1] <- n
                          y.plot.hetero.0[length(y.plot.hetero.0)+1] <- (pok*Nmax)+m 
                        }
                      }
                    }
                  }
                  coordinates <- list(x.plot.homo.1,y.plot.homo.1,x.plot.hetero.1,y.plot.hetero.1,x.plot.homo.0,y.plot.homo.0,x.plot.hetero.0,y.plot.hetero.0)
                  coordinates.length <- c()
                  for(cor in 1:length(coordinates)){
                    coordinates.length[cor] <- length(coordinates[[cor]]) 
                  }
                  coordinates.length <- max(coordinates.length)
                  for(cor in 1:length(coordinates)){
                    roznica <- coordinates.length-length(coordinates[[cor]])
                    coordinates[[cor]] <- c(coordinates[[cor]],rep(NA,times=roznica)) 
                  }
                  coordinates <- do.call(rbind,coordinates)
                  write.table(coordinates, file=file.name3,sep="\t",quote=T)
                  
                  populacja <- cbind(populacja,wektor.plci)
                  if(nrow(populacja)>1){
                    populacja <- populacja[samplex(c(1:nrow(populacja)),rep=F),]
                  }
                } else {
                  populacja <- matrix(nrow = 0,ncol = (G+1))
                }
                populacja <- rbind(populacja,selection.coeff,dominance.coeff,deparse.level = 0)
                write.table(populacja, file=file.name, sep="\t", quote=T)
              }
              if(pok == l.pokolen & rep != 1){
                file.remove(file.name3)
                file.remove(file.name4)
              }
              wrzutka
            }
            stopCluster(cl)
            
            parametry.pokolenia<-do.call(rbind,parametry.pokolenia)
            if(length(parametry.pokolenia)>0){
              colnames(parametry.pokolenia) <- c("mean_viability","heterozygosity","l.alleli","heterozygous_mutations","homozygous_mutations","fixed_loci","segregating_loci","mutations_per_individual","lethal_eq","alive_offspring")
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
  }
}

