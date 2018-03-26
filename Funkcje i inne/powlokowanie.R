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

# rozmiar.szalki <- 4
# podzial.szalki <- 2
# ruch <- 10
# l.samic <- 10
# l.samcow <- 8

powlokowanie <- function(rozmiar.szalki,podzial.szalki,ruch = 5,l.samic = 5,l.samcow = 5){
  mozliwe.ruchy <- function(x,y,z,podz.szalki=podzial.szalki,roz.szalki=rozmiar.szalki){
    ruchy <- list()
    licznik <- 0
    if(x == 1 & roz.szalki > x){
      for(k in 1:podz.szalki){
        licznik <- licznik+1
        ruchy[[licznik]] <- c(2,k)
      }
    }
    if(x == 2){
      mozliwosci <- 1
      if(podz.szalki == 1){
        #2 mozliwosci
        ruchy[[1]] <- c(1,1)
        if(roz.szalki > x){
          ruchy[[2]] <- c(3,1)
        }
      }
      if(podz.szalki == 2){
        ruchy[[1]] <- c(1,1)
        if(y == 1){
          ruchy[[2]] <- c(2,2)
        } else {ruchy[[2]] <- c(2,1)}
        licznik <- 2
        if(roz.szalki > x){
          for(k in 1:4){
            licznik <- licznik+1
            ruchy[[licznik]] <- c(3,k)
          }
        }
      }
      if(podz.szalki>2){
        mozliwosci <- mozliwosci+2+podz.szalki+2
        if(y == 1){
          sasiedzi <- c(podz.szalki,2)
          pietro <- c(podz.szalki^x,1:(podz.szalki+1))
        } else if(y == podz.szalki){
          sasiedzi <- c(y-1,1)
          pietro <- c(((podz.szalki^x)-podz.szalki):(podz.szalki^x),1)
        } else {sasiedzi <- c(y-1,y+1)
        pietro <- c(((y*podz.szalki)-podz.szalki):((y*podz.szalki)+1))}
        ruchy[[1]] <- c(1,1)
        licznik <- 1
        for(k in 1:(mozliwosci-1)){
          licznik <- licznik+1
          if(k == 1 | k == 2){
            ruchy[[licznik]] <- c(2,sasiedzi[k]) 
          } else if(k > 2 & roz.szalki > x) {
            ruchy[[licznik]] <- c(3,pietro[k-2])
          }
        }
      }
    }
    if(x > 2){
      if(podz.szalki == 1){
        #2 mozliwosci
        ruchy[[1]] <- c(x-1,1)
        if(roz.szalki > x){
          ruchy[[2]] <- c(x+1,1)
        }
      } else {
        fragmenty <- (podz.szalki)^(x-1)
        if(y %% podz.szalki > 0){
          licznik <- licznik+1
          ruchy[[licznik]] <- c(x-1,(y %/% podz.szalki)+1)
        } else if(y %% podz.szalki == 0){
          licznik <- licznik+2
          ruchy[[licznik-1]] <- c(x-1,y %/% podz.szalki)
          if(y == fragmenty){
            ruchy[[licznik]] <- c(x-1,1)
          } else {ruchy[[licznik]] <- c(x-1,(y %/% podz.szalki)+1)}
        }
        if(y %% podz.szalki == 1){
          licznik <- licznik+1
          if(y == 1){
            ruchy[[licznik]] <- c(x-1,fragmenty/podz.szalki)
          } else {ruchy[[licznik]] <- c(x-1,(y %/% podz.szalki))}
        }
        if(y == 1){
          sasiedzi <- c(fragmenty,2)
          pietro <- c(podz.szalki^x,1:(podz.szalki+1))
        } else if(y == fragmenty){
          sasiedzi <- c(y-1,1)
          pietro <- c(((podz.szalki^x)-podz.szalki):(podz.szalki^x),1)
        } else {sasiedzi <- c(y-1,y+1)
        pietro <- c(((y*podz.szalki)-podz.szalki):((y*podz.szalki)+1))}
      }
      licznik <- licznik+2
      ruchy[[licznik-1]] <- c(x,sasiedzi[1])
      ruchy[[licznik]] <- c(x,sasiedzi[2])
      if(roz.szalki>x){
        for(k in 1:(podz.szalki+2)){
          licznik <- licznik+1
          ruchy[[licznik]] <- c(x+1,pietro[k])
        } 
      }
    }
    ruchy <- do.call(rbind,ruchy)
    ruchy <- cbind(ruchy,z)
    ruchy
  }
  odleglosc <- function(samiec,samica){
    wrzutka <- 0
    pierwszy.krok <- matrix(nrow = 1,ncol = 3)
    wskaznik <- 1
    while.check <- F
    if(identical(samiec,samica)==F){
      wskaznik <- c()
      macierze <- list()
      samiec <- matrix(c(samiec,NA),nrow = 1,ncol = 3)
      while.start <- T
      while(while.start == T){
        for(i in 1:nrow(samiec)){
          macierze[[i]] <- mozliwe.ruchy(x=samiec[i,1],y=samiec[i,2],z=samiec[i,3])
        }
        samiec <- do.call(rbind,macierze)
        if(wrzutka == 0){
          samiec[,3] <- c(1:nrow(samiec))
          pierwszy.krok <- samiec
        }
        wrzutka <- wrzutka+1
        for(i in 1:nrow(samiec)){
          if(sum(as.vector(samiec[i,1:2])==samica) == 2){
            wskaznik <- c(wskaznik,as.vector(samiec[i,3]))
            while.check <- T
          }
        }
        if(while.check == T){
          wskaznik <- samplex(wskaznik,size=1)
          while.start <- F
        }
      }
    }
    wyrzut <- as.vector(c(wrzutka,pierwszy.krok[wskaznik,1:2]))
    wyrzut
  }
  if(rozmiar.szalki <= 0 | rozmiar.szalki %% 1 != 0){
    stop("rozmiar szalki musi byc liczba naturalna wieksza od 0")
  }
  if(podzial.szalki <= 0 | podzial.szalki %% 1 != 0){
    stop("rozmiar szalki musi byc liczba naturalna wieksza od 0")
  }
  if(l.samic <= 0 | l.samic %% 1 != 0 | l.samcow <= 0 | l.samcow %% 1 != 0){
    stop("liczby osobnikow musza byc liczbami naturalnymi wiekszymi od 0")
  }
  l.powlok <- 1
  if(rozmiar.szalki > 1){
    for(i in 1:(rozmiar.szalki-1)){
      if(i == 1){
        l.powlok <- c(l.powlok,podzial.szalki)
      } else {
        l.powlok <- c(l.powlok,l.powlok[i]*podzial.szalki)
      }
    }
  }
  kawalki <- matrix(nrow = sum(l.powlok),ncol=2)
  a <- 0
  b <- 1
  for(i in 1:length(l.powlok)){
    a <- a+l.powlok[i]
    kawalki[b:a,1] <- i
    kawalki[b:a,2] <- c(1:l.powlok[i])
    b <- a+1
  }
  miejsca.samic <- list()
  czarna.lista2 <- list()
  for(i in 1:l.samic){
    nr <- samplex(1:nrow(kawalki),size = 1)
    miejsca.samic[[i]] <- c(i,kawalki[nr,])
    czarna.lista2[[i]] <- NA
  }
  miejsca.samcow <- list()
  czarna.lista1 <- list()
  for(i in 1:l.samcow){
    nr <- samplex(1:nrow(kawalki),size = 1)
    miejsca.samcow[[i]] <- c(i,kawalki[nr,])
    czarna.lista1[[i]] <- NA
  }
  zakres.ruchu.samice <- rpois(l.samic,lambda = ruch)
  zakres.ruchu.samce <- rpois(l.samcow,lambda = ruch)
  kopulacje <- matrix(nrow=0,ncol=2)
  colnames(kopulacje) <- c("samiec","samica")
  if(max(c(zakres.ruchu.samice,zakres.ruchu.samce))>=1){
    for(i in 1:max(c(zakres.ruchu.samice,zakres.ruchu.samce))){
      miejsca.samic <- samplex(miejsca.samic)
      miejsca.samcow <- samplex(miejsca.samcow)
      kolejnosc <- samplex(c(0,1),size = 1)
      for(j in 1:max(c(l.samic,l.samcow))){
        if(kolejnosc == 1){
          if(l.samcow >= j){
            if(i <= zakres.ruchu.samce[miejsca.samcow[[j]][1]]){
              nr.samca <- miejsca.samcow[[j]][1]
              samiec.place <- miejsca.samcow[[j]][2:3]
              wybranki <- c(1:l.samic)
              wyb <- do.call(rbind,miejsca.samic)
              for(l in 1:length(czarna.lista1[[nr.samca]])){
                if(sum(is.na(czarna.lista1[[nr.samca]]))==0){
                  num <- which(wyb[,1]==czarna.lista1[[nr.samca]][l])
                  num1 <- which(wybranki == num)
                  wybranki <- wybranki[-num1]
                }
              }
              if(length(wybranki)>0){
                drogi <- matrix(nrow = length(wybranki),ncol = 4)
                for(k in 1:length(wybranki)){
                  drogi[k,] <- c(odleglosc(samiec.place,miejsca.samic[[wybranki[k]]][2:3]),miejsca.samic[[wybranki[k]]][1])
                }
                if(sum(drogi[,1]==0) != 0){
                  nr <- which(drogi[,1] == 0)
                  nr <- samplex(nr,size = 1)
                  kopulacje <- rbind(kopulacje,c(nr.samca,drogi[nr,4]))
                  czarna.lista1[[nr.samca]] <- as.vector(na.omit(c(czarna.lista1[[nr.samca]],drogi[nr,4])))
                  czarna.lista2[[drogi[nr,4]]] <- as.vector(na.omit(c(czarna.lista2[[drogi[nr,4]]],nr.samca)))
                } else {
                  drogi <- drogi[samplex(1:nrow(drogi)),]
                  if(is.null(nrow(drogi))==T){
                    drogi <- matrix(drogi,nrow = 1,ncol = 4)
                  }
                  nr <- which.min(drogi[,1])
                  miejsca.samcow[[j]][2:3] <- drogi[nr,2:3]
                }
              }
            }
          }

          if(l.samic >= j){
            if(i <= zakres.ruchu.samice[miejsca.samic[[j]][1]]){
              nr.samicy <- miejsca.samic[[j]][1]
              samica.place <- miejsca.samic[[j]][2:3]
              wyb <- do.call(rbind,miejsca.samcow)
              wybranki <- c(1:l.samcow)
              for(l in 1:length(czarna.lista2[[nr.samicy]])){
                if(sum(is.na(czarna.lista2[[nr.samicy]]))==0){
                  num <- which(wyb[,1]==czarna.lista2[[nr.samicy]][l])
                  num1 <- which(wybranki == num)
                  wybranki <- wybranki[-num1]
                }
              }
              if(length(wybranki)>0){
                drogi <- matrix(nrow = length(wybranki),ncol = 4)
                for(k in 1:length(wybranki)){
                  drogi[k,] <- c(odleglosc(samica.place,miejsca.samcow[[wybranki[k]]][2:3]),miejsca.samcow[[wybranki[k]]][1])
                }
                if(sum(drogi[,1]==0) != 0){
                  nr <- which(drogi[,1] == 0)
                  nr <- samplex(nr,size = 1)
                  kopulacje <- rbind(kopulacje,c(drogi[nr,4],nr.samicy))
                  czarna.lista2[[nr.samicy]] <- as.vector(na.omit(c(czarna.lista2[[nr.samicy]],drogi[nr,4])))
                  czarna.lista1[[drogi[nr,4]]] <- as.vector(na.omit(c(czarna.lista1[[drogi[nr,4]]],nr.samicy)))
                } else {
                  drogi <- drogi[samplex(1:nrow(drogi)),]
                  if(is.null(nrow(drogi))==T){
                    drogi <- matrix(drogi,nrow = 1,ncol = 4)
                  }
                  nr <- which.min(drogi[,1])
                  miejsca.samic[[j]][2:3] <- drogi[nr,2:3]
                }
              }
            }
          }
        } else {
          if(l.samic >= j){
            if(i <= zakres.ruchu.samice[miejsca.samic[[j]][1]]){
              nr.samicy <- miejsca.samic[[j]][1]
              samica.place <- miejsca.samic[[j]][2:3]
              wyb <- do.call(rbind,miejsca.samcow)
              wybranki <- c(1:l.samcow)
              for(l in 1:length(czarna.lista2[[nr.samicy]])){
                if(sum(is.na(czarna.lista2[[nr.samicy]]))==0){
                  num <- which(wyb[,1]==czarna.lista2[[nr.samicy]][l])
                  num1 <- which(wybranki == num)
                  wybranki <- wybranki[-num1]
                }
              }
              if(length(wybranki)>0){
                drogi <- matrix(nrow = length(wybranki),ncol = 4)
                for(k in 1:length(wybranki)){
                  drogi[k,] <- c(odleglosc(samica.place,miejsca.samcow[[wybranki[k]]][2:3]),miejsca.samcow[[wybranki[k]]][1])
                }
                if(sum(drogi[,1]==0) != 0){
                  nr <- which(drogi[,1] == 0)
                  nr <- samplex(nr,size = 1)
                  kopulacje <- rbind(kopulacje,c(drogi[nr,4],nr.samicy))
                  czarna.lista2[[nr.samicy]] <- as.vector(na.omit(c(czarna.lista2[[nr.samicy]],drogi[nr,4])))
                  czarna.lista1[[drogi[nr,4]]] <- as.vector(na.omit(c(czarna.lista1[[drogi[nr,4]]],nr.samicy)))
                } else {
                  drogi <- drogi[samplex(1:nrow(drogi)),]
                  if(is.null(nrow(drogi))==T){
                    drogi <- matrix(drogi,nrow = 1,ncol = 4)
                  }
                  nr <- which.min(drogi[,1])
                  miejsca.samic[[j]][2:3] <- drogi[nr,2:3]
                }
              }
            }
          }
          if(l.samcow >= j){
            if(i <= zakres.ruchu.samce[miejsca.samcow[[j]][1]]){
              nr.samca <- miejsca.samcow[[j]][1]
              samiec.place <- miejsca.samcow[[j]][2:3]
              wybranki <- c(1:l.samic)
              wyb <- do.call(rbind,miejsca.samic)
              for(l in 1:length(czarna.lista1[[nr.samca]])){
                if(sum(is.na(czarna.lista1[[nr.samca]]))==0){
                  num <- which(wyb[,1]==czarna.lista1[[nr.samca]][l])
                  num1 <- which(wybranki == num)
                  wybranki <- wybranki[-num1]
                }
              }
              if(length(wybranki)>0){
                drogi <- matrix(nrow = length(wybranki),ncol = 4)
                for(k in 1:length(wybranki)){
                  drogi[k,] <- c(odleglosc(samiec.place,miejsca.samic[[wybranki[k]]][2:3]),miejsca.samic[[wybranki[k]]][1])
                }
                if(sum(drogi[,1]==0) != 0){
                  nr <- which(drogi[,1] == 0)
                  nr <- samplex(nr,size = 1)
                  kopulacje <- rbind(kopulacje,c(nr.samca,drogi[nr,4]))
                  czarna.lista1[[nr.samca]] <- as.vector(na.omit(c(czarna.lista1[[nr.samca]],drogi[nr,4])))
                  czarna.lista2[[drogi[nr,4]]] <- as.vector(na.omit(c(czarna.lista2[[drogi[nr,4]]],nr.samca)))
                } else {
                  drogi <- drogi[samplex(1:nrow(drogi)),]
                  if(is.null(nrow(drogi))==T){
                    drogi <- matrix(drogi,nrow = 1,ncol = 4)
                  }
                  nr <- which.min(drogi[,1])
                  miejsca.samcow[[j]][2:3] <- drogi[nr,2:3]
                }
              }
            }
          }
        }
      }
    }
  }
  kopulacje
}
