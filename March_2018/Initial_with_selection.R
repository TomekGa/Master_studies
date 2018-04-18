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
    if(length(l.chromosomes)>1){
      for(a in 1:(length(l.chromosomes)-1)){
        boundaries[a] <- l.chromosomes[a] + boundaries.temporary[a]
        boundaries.temporary[a+1] <- boundaries[a]
      }
      # number of chiasms on each chromosome
      n.chiasms <- rpois((length(boundaries)+1),lambda = lambda)
    } else {
      boundaries <- NULL
      n.chiasms <- rpois(1,lambda = lambda)
    }
  }

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

# Next step----------------------------------------------------------------

N<-1000 #size of populations (even number)
M<-2900 #size of haploid genome uder selection (in loci, 2 chromosomes)
L<-2 #average number of chiasms on each chromosome
G<-3000 #overall size of haploid genome (in loci, 200 neutral loci included)
#Ul<-0.03 #lethal mutations rate (from Wang's article)
U<-0.5 #deleterious mutations rate (Simmons, 1977)
l.pokolen <- 2
l.replikatow <- 1
l.plemnikow<-100
l.jaj<-10
s<-0.05 #mean selection coefficient (Halligan, Keightley, 2009)
h<-0.36 #mean dominance coefficient (Garcia-Dorado,Caballero, 2001)
u<-U/(2*M) #average rate of mutation per locus
neutral.allels.nr <- 7
choice.coeff.wektor <- 1.3
no_cores <- 3

Nmax <- N
N.fem<-N/2
N.mal<-N/2
loci.neutralne<-c(1:100)*30
wszystkie.loci<-c(1:G)
loci.sel<-wszystkie.loci[!is.element(wszystkie.loci, loci.neutralne)]
LETTERSplus <- LETTERS[-which(LETTERS=="M" | LETTERS == "N")]
allele.neutralne<-LETTERSplus[1:neutral.allels.nr]

selection.coeff <- rep(s, times = G)
selection.coeff[loci.neutralne] <- NA
selection.coeff <- c(selection.coeff,NA)

dominance.coeff <- selection.coeff
dominance.coeff[!is.na(dominance.coeff)] <- h

populacja<-matrix(nrow=Nmax, ncol=G+1)
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

library(foreach) #introduction to foreach
library(doParallel)
cl <- makeCluster(no_cores)
registerDoParallel(cl)
init.pop <- c(1:N)
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

populacja[,G+1] <- samplex(c(0,1),size = Nmax,replace = T)
populacja <- rbind(populacja,selection.coeff,dominance.coeff,deparse.level = 0)
write.table(populacja, file="populacja.txt", sep="\t", quote=T)

writeLines("","log.txt")
for(pok in 1:l.pokolen){
  start <- Sys.time()
  sink("log.txt",append = T)
  cat(paste("pokolenie",pok,"\n",sep = " "))
  sink()
  
  populacja <- read.table("populacja.txt",sep = "\t",header = T,colClasses = "character")
  selection.coeff <- as.numeric(populacja[(nrow(populacja)-1),])
  dominance.coeff <- as.numeric(populacja[nrow(populacja),])
  populacja <- populacja[1:(nrow(populacja)-2),]
  N <- nrow(populacja)
  samce <- samplex(which(populacja[,(G+1)]==1))
  samice <- samplex(c(1:N)[!is.element(c(1:N),samce)],replace = F)
  populacja <- populacja[,(1:G)]
  
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
  
  chromosomy <- matrix(nrow=2*N, ncol=G) #2 rows per one individual
  for(i in 1:N){
    for(j in 1:G){
      chromosomy[2*i-1,j]<-substring(populacja[i,j],0,1)
      chromosomy[2*i,j]<-substring(populacja[i,j],2)
    }
  }
  
  plemniki<-matrix(nrow=(length(samce)*l.plemnikow),ncol=G)
  jaja<-matrix(nrow=(length(samice)*l.jaj),ncol=G)
  
  for(i in 1:length(samice)){
    nic2<-samice[i]
    C <- chromosomy[(2*nic2-1),]
    D <- chromosomy[(2*nic2),]
    for(nr in 1:l.jaj){
      jaja[((i-1)*l.jaj)+nr,]<-crossing.over(C,D,l.chromosomes = G)
      
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
      plemniki[((i-1)*l.plemnikow)+nr,]<-crossing.over(A,B,l.chromosomes = G)
      ranking.plemnikow[((i-1)*l.plemnikow)+nr] <- sum(plemniki[((i-1)*l.plemnikow)+nr,] == "M")
    }
  }
  
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
    miejsce.plemnika <- samplex(c(1:nrow(plemniki)),size=1,prob = ranking.prob)
    gameta1 <- plemniki[miejsce.plemnika,1:G]
    plemniki <- plemniki[-miejsce.plemnika,]
    for(x in 1:G){
      populacja1[y,x] <- paste(gameta1[x],gameta2[x],sep="")
    }
  }
  
  populacja1 <- cbind(populacja1,c(samplex(c(0,1),size=nrow(populacja1),replace = T)))
  
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
  wektor.plci <- new_population[,(G+2)]
  populacja<-new_population[,2:(G+1)]
  if(is.null(nrow(populacja))==T){
    populacja <- matrix(populacja,nrow = 1,ncol=(length(populacja)))
  }
  
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

  populacja <- cbind(populacja,wektor.plci)
  populacja <- rbind(populacja,selection.coeff,dominance.coeff,deparse.level = 0)
  write.table(populacja, file="populacja.txt", sep="\t", quote=T)
  print(Sys.time()-start)
}















