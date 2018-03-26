zmutowana.gameta <- function(x, wspolczynnik){
  if(x %% 1 != 0 | x < 0){
    stop("Liczba mutacji musi byc liczba NATURALNA wieksza badz rowna 0")
  }
  if(wspolczynnik <= 0){
    stop("wspolczynnik musi byc wiekszy od zera")
  }
  if(wspolczynnik < 1){
    warning("odwrotna zaleznosc - skrajne liczby mutacji sa bardziej prawdopodobne")
  }
  if(wspolczynnik == 1){
    warning("prawdopodobienstwo kazdej liczby mutacji jest rowne")
  }
  ### wirtualny histogram
  l.slupkow <- x+1
  l.cegiel.slupek <- c()
  numery.slupkow <- c(0:(l.slupkow-1))
  if(l.slupkow == 1){
    l.cegiel.slupek[1] <- 1
  }
  if(l.slupkow > 1){
    if(l.slupkow %% 2 != 0){
      srodek <- sum(numery.slupkow)/l.slupkow
      l.cegiel.slupek[srodek+1] <- 1
      for (i in 1:srodek){
        l.cegiel.slupek[srodek+1-i] <- l.cegiel.slupek[(srodek+1-i)+1]/wspolczynnik
        l.cegiel.slupek[srodek+1+i] <- l.cegiel.slupek[(srodek+1+i)-1]/wspolczynnik
      }
    } else {
      srodek <- sum(c(0:numery.slupkow[l.slupkow-1]))/(l.slupkow-1)
      srodek[2] <- srodek[1]+1
      l.cegiel.slupek[srodek+1] <- 1
      l.stopni <- srodek[1]
      if(l.stopni > 0){
        for (i in 1:l.stopni){
          l.cegiel.slupek[srodek[1]+1-i] <- l.cegiel.slupek[(srodek[1]+1-i)+1]/wspolczynnik
          l.cegiel.slupek[srodek[2]+1+i] <- l.cegiel.slupek[(srodek[2]+1+i)-1]/wspolczynnik
        }
      }
    }
  }
  l.cegielek <- sum(l.cegiel.slupek)
  prawdopodobienstwo.cegielki <- 1/l.cegielek
  prawdopodobienstwa <- l.cegiel.slupek * prawdopodobienstwo.cegielki
  ###Tworzenie gamety
  gameta <- sample(c(0:x),size = 1,prob = prawdopodobienstwa)
  gameta
}
liczba.potomstwa <- function(x,wspolczynnik,max.potomstwa,prog.smiertelnosci,podkrecanie = FALSE){
  if(wspolczynnik <= 1){
    stop("argument _wspolczynnik_ musi być większy od 1")
  }
  odwrotnosc <- prog.smiertelnosci-x
  l.cegielek <- 1
  for(i in 1:prog.smiertelnosci){
    l.cegielek[i+1] <- l.cegielek[i]*wspolczynnik
  }
  l.cegielek <- l.cegielek[1:prog.smiertelnosci]
  l.cegielek <- c(0,l.cegielek)
  cegielka <- max.potomstwa/l.cegielek[prog.smiertelnosci+1]
  l.potomstwa <- round(l.cegielek[odwrotnosc+1]*cegielka)
  if(podkrecanie == T & x < prog.smiertelnosci & l.potomstwa == 0){
    l.potomstwa <- 1
  }
  l.potomstwa
}

N <- 1000
l.pokolen <- 1000
prog.smiertelnosci <- 100 #musi być parzysta
l.potomstwa.max <- 100
l.jaj.max <- 100
l.plemnikow.max <- 100
srednia.mutacji <- 1000
srednia.mutacji.pozytyw <- 1
wspolczynnik1 <- seq(1.01,1.5,by=0.01)
wspolczynnik2 <- seq(0.1,5,by = 0.1)
constants <- paste("N",N,"Generations",l.pokolen,"Number of offspring",l.potomstwa.max,"Mutation mean",srednia.mutacji,"Possitive mutation mean",srednia.mutacji.pozytyw)
#par(mfrow=c(2,2))
#par(mfrow=c(1,1))
writeLines("","log.txt")
#writeLines(c(constants,"At least one population became extinct"),"junk.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations higher than 5"),"smaller_junk.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations between 4 and 5"),"diff_5.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations between 3 and 4"),"diff_4.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations between 2 and 3"),"diff_3.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations between 1 and 2"),"diff_2.txt")
writeLines(c(constants,"Difference in mean number of mutations in last 50 generations below 1"),"diff_1.txt")
library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

invisible(foreach(iter = wspolczynnik1) %:%
  foreach(iter2 = wspolczynnik2) %dopar% {
  sink("log.txt",append = T)
  cat(paste(iter,iter2,"\n"))
  sink()
  # Populacja aseksualna ----------------------------------------------------
  populacja <- sample(0:(prog.smiertelnosci-1),size=N,replace = TRUE)
  
  srednia.ilosc.mutacji <- c()
  for(l in 1:l.pokolen){
    #print(l)
    NN <- length(populacja)
    #print(NN)
    l.mutacji <- rpois(1,srednia.mutacji)
    if(l.mutacji >= 1){
      miejsca.mutacji <- sample(1:NN,size=l.mutacji,replace = T)
      for(i in miejsca.mutacji){
        populacja[i] <- populacja[i]+1
      }
    }
    l.mutacji.pozytyw <- rpois(1,srednia.mutacji.pozytyw)
    if(l.mutacji.pozytyw >= 1){
      miejsca.mutacji.pozytyw <- sample(1:NN,size=l.mutacji.pozytyw,replace = T)
      for(i in miejsca.mutacji.pozytyw){
        populacja[i] <- populacja[i]-1
        if(populacja[i]<0){populacja[i] <- 0}
      }
    }
    populacja <- populacja[!populacja >= prog.smiertelnosci]
    if(length(populacja)==0){
      srednia.ilosc.mutacji[l:l.pokolen] <- prog.smiertelnosci
      break
    }
    l.potomstwa <- c()
    for(p in 1:length(populacja)){
      l.potomstwa[p] <- liczba.potomstwa(populacja[p],wspolczynnik = iter,max.potomstwa = l.potomstwa.max,prog.smiertelnosci = prog.smiertelnosci,podkrecanie = T)
    }
    l.potomstwa <- c(0,l.potomstwa)
    for(m in 1:length(populacja)){
      l.potomstwa[m+1] <- l.potomstwa[m] + l.potomstwa[m+1]
    }
    populacja1 <- c()
    for(k in 1:length(populacja)){
      populacja1[(l.potomstwa[k]+1):l.potomstwa[k+1]] <- populacja[k]
    }
    if(length(populacja1)>=N){
      populacja1 <- sample(populacja1,size = N,replace = F)
    }else {populacja1 <- sample(populacja1,size = length(populacja1),replace = F)}
    
    srednia.ilosc.mutacji[l] <- mean(populacja1)
    populacja <- populacja1
  }
  # plot(c(1:l.pokolen),srednia.ilosc.mutacji,main="1.3")
  # fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
  # co <- coef(fit)
  # abline(fit, col="blue", lwd=2)
  
  # Populacja seksualna z taką samą selekcją --------------------------------
  
  if(srednia.ilosc.mutacji[l.pokolen]==prog.smiertelnosci){
    srednia.ilosc.mutacji.2 <- c()
    srednia.ilosc.mutacji.2[1:l.pokolen] <- prog.smiertelnosci
  } else {
    populacja <- sample(0:(prog.smiertelnosci-1),size=N,replace = TRUE)
    wektor.plci <- sample(c(0,1),size = length(populacja),replace = T)
    
    srednia.ilosc.mutacji.2 <- c()
    for(l in 1:l.pokolen){
      #print(l)
      NN <- length(populacja)
      l.mutacji <- rpois(1,srednia.mutacji)
      if(l.mutacji >= 1){
        miejsca.mutacji <- sample(1:NN,size=l.mutacji,replace = T)
        for(i in miejsca.mutacji){
          populacja[i] <- populacja[i]+1 
        }
      }
      l.mutacji.pozytyw <- rpois(1,srednia.mutacji.pozytyw)
      if(l.mutacji.pozytyw >= 1){
        miejsca.mutacji.pozytyw <- sample(1:NN,size=l.mutacji.pozytyw,replace = T)
        for(i in miejsca.mutacji.pozytyw){
          populacja[i] <- populacja[i]-1
          if(populacja[i]<0){populacja[i] <- 0}
        }
      }
      if(sum(populacja >= prog.smiertelnosci) != 0){
        ponad.prog <- which(populacja >= prog.smiertelnosci)
        populacja <- populacja[-ponad.prog]
        wektor.plci <- wektor.plci[-ponad.prog]
      }
      if(length(populacja)==0){
        srednia.ilosc.mutacji.2[l:l.pokolen] <- prog.smiertelnosci
        break
      }
      samce <- which(wektor.plci == 0)
      samice <- which(wektor.plci == 1)
      
      l.jaj <- c()
      l.plemnikow <- c()
      for(p in 1:length(samice)){
        l.jaj[p] <- liczba.potomstwa(populacja[samice[p]],iter,max.potomstwa = l.jaj.max,prog.smiertelnosci = prog.smiertelnosci,podkrecanie = T)
      }
      for(p in 1:length(samce)){
        l.plemnikow[p] <- liczba.potomstwa(populacja[samce[p]],iter,max.potomstwa = l.plemnikow.max,prog.smiertelnosci = prog.smiertelnosci, podkrecanie = T)
      }
      l.jaj <- c(0,l.jaj)
      for(m in 1:length(samice)){
        l.jaj[m+1] <- l.jaj[m] + l.jaj[m+1]
      }
      l.plemnikow <- c(0,l.plemnikow)
      for(m in 1:length(samce)){
        l.plemnikow[m+1] <- l.plemnikow[m] + l.plemnikow[m+1]
      }
      jaja <- c()
      for(k in 1:length(samice)){
        jaja[(l.jaj[k]+1):l.jaj[k+1]] <- zmutowana.gameta(populacja[samice[k]],wspolczynnik = iter2)
      }
      jaja <- sample(jaja)
      #print(length(jaja))
      plemniki <- c()
      for(k in 1:length(samce)){
        plemniki[(l.plemnikow[k]+1):l.plemnikow[k+1]] <- zmutowana.gameta(populacja[samce[k]],wspolczynnik = iter2)
      }
      plemniki <- sample(plemniki)
      #print(length(plemniki))
      populacja1 <- c()
      if(length(jaja)<=length(plemniki)){
        plemniki <- sample(plemniki,size = length(jaja),replace = F)
        populacja1 <- jaja + plemniki
      }
      else { 
        jaja <- sample(jaja,size = length(plemniki),replace = F)
        populacja1 <- jaja + plemniki
      }
      wektor.plci <- sample(c(0,1),size=length(populacja1),replace = T)
      
      if(length(populacja1)>=N){
        wybrancy <- sample(c(1:length(populacja1)),size = N,replace = F)
        populacja1 <- populacja1[wybrancy]
        wektor.plci <- wektor.plci[wybrancy]
      }
      srednia.ilosc.mutacji.2[l] <- mean(populacja1)
      populacja <- populacja1
    }
    
    # plot(c(1:l.pokolen),srednia.ilosc.mutacji)
    # fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
    # co <- coef(fit)
    # abline(fit, col="blue", lwd=2)
  }
# wrzucanie do plików -----------------------------------------------------
  smut1 <- paste(srednia.ilosc.mutacji,collapse = " ")
  smut2 <- paste(srednia.ilosc.mutacji.2,collapse = " ")
  wrzutka <- paste("Coef. # of offspring",iter,"Coef. linked mutation",iter2,"\n","asex",smut1,"\n","sex",smut2,"\n")
  roznica <- abs(mean(srednia.ilosc.mutacji[(l.pokolen-49):l.pokolen])-mean(srednia.ilosc.mutacji.2[(l.pokolen-49):l.pokolen]))
  if(mean(srednia.ilosc.mutacji[(l.pokolen-49):l.pokolen])==prog.smiertelnosci | mean(srednia.ilosc.mutacji.2[(l.pokolen-49):l.pokolen])==prog.smiertelnosci){
    # sink("junk.txt",append = T)
    # cat(wrzutka)
    # sink()
  } else if(roznica>5){
    sink("smaller_junk.txt",append = T)
    cat(wrzutka)
    sink()
  } else if(roznica>4){
    sink("diff_5.txt",append = T)
    cat(wrzutka)
    sink()
  } else if(roznica>3){
    sink("diff_4.txt",append = T)
    cat(wrzutka)
    sink()
  } else if(roznica>2){
    sink("diff_3.txt",append = T)
    cat(wrzutka)
    sink()
  } else if(roznica>1){
    sink("diff_2.txt",append = T)
    cat(wrzutka)
    sink()
  } else {
    sink("diff_1.txt",append = T)
    cat(wrzutka)
    sink()
  }
})
stopCluster(cl)
print("Finished")



