zmutowana.gameta <- function(x, wspolczynnik){
  if(x %% 1 != 0 | x < 0){
    stop("Liczba mutacji musi byæ liczb¹ NATURALN¥ wiêksz¹ b¹dŸ równ¹ 0")
  }
  if(wspolczynnik <= 0){
    stop("wspolczynnik musi byæ wiêkszy od zera")
  }
  if(wspolczynnik < 1){
    warning("odwrotna zale¿noœæ - skrajne liczby mutacji s¹ bardziej prawdopodobne")
  }
  if(wspolczynnik == 1){
    warning("prawdopodobieñstwo ka¿dej liczby mutacji jest równe")
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
    stop("argument _wspolczynnik_ musi byæ wiêkszy od 1")
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

N <- 10000
l.pokolen <- 1000
prog.smiertelnosci <- 100 #musi byæ parzysta
l.potomstwa.max <- 50
l.jaj.max <- 500
l.plemnikow.max <- 500
srednia.mutacji <- 1000
srednia.mutacji.pozytyw <- 1
wspolczynnik1 <- seq(1.01,1.5,by=0.01)
wspolczynnik2 <- seq(0.1,5,by = 0.1)
#par(mfrow=c(2,2))
par(mfrow=c(1,1))
# Populacja aseksualna ----------------------------------------------------
populacja <- sample(0:(prog.smiertelnosci-1),size=N,replace = TRUE)

srednia.ilosc.mutacji <- c()
for(l in 1:l.pokolen){
  print(l)
  NN <- length(populacja)
  print(NN)
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
    l.potomstwa[p] <- liczba.potomstwa(populacja[p],wspolczynnik = 1.3,max.potomstwa = l.potomstwa.max,prog.smiertelnosci = prog.smiertelnosci,podkrecanie = T)
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
plot(c(1:l.pokolen),srednia.ilosc.mutacji,main="1.3")
fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
co <- coef(fit)
abline(fit, col="blue", lwd=2)

# Populacja seksualna z tak¹ sam¹ selekcj¹ --------------------------------

populacja <- sample(0:(prog.smiertelnosci-1),size=N,replace = TRUE)
wektor.plci <- sample(c(0,1),size = length(populacja),replace = T)
  
srednia.ilosc.mutacji <- c()
for(l in 1:l.pokolen){
  print(l)
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
    srednia.ilosc.mutacji[l:l.pokolen] <- prog.smiertelnosci
    break
  }
  samce <- which(wektor.plci == 0)
  samice <- which(wektor.plci == 1)
  
  l.jaj <- c()
  l.plemnikow <- c()
  for(p in 1:length(samice)){
    l.jaj[p] <- liczba.potomstwa(populacja[samice[p]],1.25,max.potomstwa = l.jaj.max,prog.smiertelnosci = prog.smiertelnosci)
  }
  for(p in 1:length(samce)){
    l.plemnikow[p] <- liczba.potomstwa(populacja[samce[p]],1.25,max.potomstwa = l.plemnikow.max,prog.smiertelnosci = prog.smiertelnosci)
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
    jaja[(l.jaj[k]+1):l.jaj[k+1]] <- zmutowana.gameta(populacja[samice[k]],wspolczynnik = 5.3)
  }
  jaja <- sample(jaja)
  print(length(jaja))
  plemniki <- c()
  for(k in 1:length(samce)){
    plemniki[(l.plemnikow[k]+1):l.plemnikow[k+1]] <- zmutowana.gameta(populacja[samce[k]],wspolczynnik = 5.3)
  }
  print(length(plemniki))
  populacja1 <- c()
  for(n in 1:length(jaja)){
    if(length(plemniki)>0){
      wybrany.plemnik <- sample(c(1:length(plemniki)),size=1)
      populacja1[n] <- sum(jaja[n],plemniki[wybrany.plemnik])
      plemniki <- plemniki[-wybrany.plemnik]
    }
  }
  wektor.plci <- sample(c(0,1),size=length(populacja1),replace = T)
  
  if(length(populacja1)>=N){
    wybrancy <- sample(c(1:length(populacja1)),size = N,replace = F)
    populacja1 <- populacja1[wybrancy]
    wektor.plci <- wektor.plci[wybrancy]
  }
  srednia.ilosc.mutacji[l] <- mean(populacja1)
  populacja <- populacja1
}

plot(c(1:l.pokolen),srednia.ilosc.mutacji)
fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
co <- coef(fit)
abline(fit, col="blue", lwd=2)


