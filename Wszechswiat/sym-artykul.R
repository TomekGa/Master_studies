zmutowane.dzieci <- function(x,y){
  ### wirtualny histogram
  l.slupkow <- sum(x,y)+1
  l.cegiel.slupek <- c()
  l.cegiel.slupek[1] <- 1
  numery.slupkow <- c(0:(l.slupkow-1))
  if(l.slupkow > 1){
    for(j in 2:l.slupkow){
      if((mean(c(x,y))-numery.slupkow[j-1]) > 0.5){
        l.cegiel.slupek[j] <- l.cegiel.slupek[j-1]+1
      }
      else if((mean(c(x,y))-numery.slupkow[j-1]) == 0.5){
        l.cegiel.slupek[j] <- l.cegiel.slupek[j-1]
      }
      else{l.cegiel.slupek[j] <- l.cegiel.slupek[j-1]-1}
    }
  }
  l.cegielek <- sum(l.cegiel.slupek)
  prawdopodobienstwo.cegielki <- 1/l.cegielek
  prawdopodobienstwa <- l.cegiel.slupek * prawdopodobienstwo.cegielki
  ###Tworzenie dziecka
  dziecko <- sample(c(0:sum(x,y)),size = 1,prob = prawdopodobienstwa)
  dziecko
}

### Za³o¿enia: po³owa osobników dostaje jedn¹ mutacjê oraz 1 samica ma zawsze 
# 2 potomstwa w stosunku 1:1

# Symulacja - Tylko dobór naturalny ---------------------------------------
N <- 1000
l.potomstwa <- 2
l.pokolen <- 1000
prog.smiertelnosci <- 20 #musi byæ parzysta
ancestralna <- sample(1:(prog.smiertelnosci-1),size=N,replace = TRUE)
populacja <- ancestralna
l.mutacji <- 500

srednia.ilosc.mutacji <- c()
for(l in 1:l.pokolen){
  miejsca.mutacji <- sample(1:N,size=l.mutacji,replace = F)
  for(i in miejsca.mutacji){
    populacja[i] <- populacja[i]+1 
  }
  bez.dziesiatek <- populacja[!populacja >= prog.smiertelnosci]
  for(j in 1:N){
    if(populacja[j] >= prog.smiertelnosci){
      populacja[j] <- sample(bez.dziesiatek,size=1)
    }
  }
  samce <- sample(1:N,size = (N/2), replace = F)
  samice <- sample(c(1:N)[!is.element(1:N,samce)],size = (N/2), replace = F)
  
  populacja1 <- c()
  for(k in 1:(N/2)){
    populacja1[k] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T)])
    populacja1[k+(N/2)] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T)])
  }
  srednia.ilosc.mutacji[l] <- mean(populacja1)
  populacja <- populacja1
}

plot(c(1:l.pokolen),srednia.ilosc.mutacji)
fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
co <- coef(fit)
abline(fit, col="blue", lwd=2)

# Symulacja - Dobór p³ciowy - równa szkodliwoœæ mutacji -------------------
N <- 1000
l.potomstwa <- 2
l.pokolen <- 1000
prog.smiertelnosci <- 20 #musi byæ parzysta
ancestralna <- sample(1:(prog.smiertelnosci-1),size=N,replace = TRUE)
populacja <- ancestralna
l.mutacji <- 500

srednia.ilosc.mutacji <- c()
for(l in 1:l.pokolen){
  miejsca.mutacji <- sample(1:N,size=l.mutacji,replace = F)
  for(i in miejsca.mutacji){
    populacja[i] <- populacja[i]+1 
  }
  bez.dziesiatek <- populacja[!populacja >= prog.smiertelnosci]
  for(j in 1:N){
    if(populacja[j] >= prog.smiertelnosci){
      populacja[j] <- sample(bez.dziesiatek,size=1)
    }
  }
  samce <- sample(1:N,size = (N/2), replace = F)
  samice <- sample(c(1:N)[!is.element(1:N,samce)],size = (N/2), replace = F)
  mutacje <- populacja[samce]
  for(n in 1:length(mutacje)){
    mutacje[n] <- (prog.smiertelnosci-1)-mutacje[n]
  }
  suma.mutacji <- sum(mutacje)
  cegielka <- 1/suma.mutacji
  probabilities <- samce
  for(m in 1:length(probabilities)){
    probabilities[m] <- mutacje[m]*cegielka
  }
  
  populacja1 <- c()
  for(k in 1:(N/2)){
    populacja1[k] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T,prob=probabilities)])
    populacja1[k+(N/2)] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T,prob=probabilities)])
  }
  srednia.ilosc.mutacji[l] <- mean(populacja1)
  populacja <- populacja1
}

plot(c(1:l.pokolen),srednia.ilosc.mutacji)
fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
co <- coef(fit)
abline(fit, col="blue", lwd=2)

# Dobór p³ciowy - 2x bardziej szkodliwe na samce --------------------------
N <- 1000
l.potomstwa <- 2
l.pokolen <- 1000
prog.smiertelnosci <- 20 #musi byæ parzysta
prog.samcow <- prog.smiertelnosci/2
ancestralna <- sample(1:(prog.smiertelnosci-1),size=N,replace = TRUE)
populacja <- ancestralna
l.mutacji <- 500

srednia.ilosc.mutacji <- c()
for(l in 1:l.pokolen){
  samce <- sample(1:N,size = (N/2), replace = F)
  samice <- sample(c(1:N)[!is.element(1:N,samce)],size = (N/2), replace = F)
  miejsca.mutacji <- sample(1:N,size=l.mutacji,replace = F)
  for(i in miejsca.mutacji){
    populacja[i] <- populacja[i]+1 
  }
  bez.samcow <- populacja[samce] < prog.samcow
  chore.samce <- samce[!bez.samcow]
  zdrowe.samce <- !is.element(samce,chore.samce)
  samce <- samce[zdrowe.samce]
  zdrowe.os <- !is.element(c(1:N),chore.samce)
  bez.dziesiatek <- populacja[populacja < prog.smiertelnosci & zdrowe.os]
  for(j in 1:N){
    if(populacja[j] >= prog.smiertelnosci | zdrowe.os[j] == F){
      populacja[j] <- sample(bez.dziesiatek,size=1)
    }
  }

  mutacje <- populacja[samce]
  for(n in 1:length(mutacje)){
    mutacje[n] <- (prog.smiertelnosci-1)-mutacje[n]
  }
  suma.mutacji <- sum(mutacje)
  cegielka <- 1/suma.mutacji
  probabilities <- samce
  for(m in 1:length(probabilities)){
    probabilities[m] <- mutacje[m]*cegielka
  }
  
  populacja1 <- c()
  for(k in 1:(N/2)){
    populacja1[k] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T,prob=probabilities)])
    populacja1[k+(N/2)] <- zmutowane.dzieci(populacja[samice[k]],populacja[sample(samce,size = 1,replace = T,prob=probabilities)])
  }
  srednia.ilosc.mutacji[l] <- mean(populacja1)
  populacja <- populacja1
}

plot(c(1:l.pokolen),srednia.ilosc.mutacji)
fit <- glm(srednia.ilosc.mutacji~c(1:l.pokolen))
co <- coef(fit)
abline(fit, col="blue", lwd=2)

