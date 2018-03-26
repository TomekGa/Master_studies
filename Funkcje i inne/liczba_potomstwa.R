
# Liczba potomstwa --------------------------------------------------------

#wpolczynnik - ile razy wiêcej potomstwa ma osobnik z 1 mutacj¹ mniej
#podkrecanie - dla liczby mutacji poni¿ej progu œmiertelnoœci liczba potomstwa nie mo¿e spaœæ poni¿ej 1
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

# skrypt do rysowania wykresu liczby potomstwa
wspolczynniki <- c(seq(1.01,1.1,by=0.01),1.15,1.2,1.3,1.5)
iksy <- c(0:100)
y <- c()
for(j in wspolczynniki){
  igreki <- c()
  for(k in 0:100){
    igreki[k+1] <- liczba.potomstwa(k,wspolczynnik = j,max.potomstwa = 100,prog.smiertelnosci = 100,podkrecanie = T) 
  }
  if(j == wspolczynniki[1]){
    plot(iksy,igreki,type="l",yaxp = c(0,100,50),ylab = "liczba potomstwa",xlab = "liczba mutacji",xaxp = c(0,100,50),main = "Liczba potomstwa w zale¿noœci od wspó³czynnika \n max potomstwa = 100, próg bezp³odnoœci = 100")
    
  } else {
    lines(iksy,igreki)
  }
  y[which(j == wspolczynniki)] <- igreki[11]+1.5
}
x <- rep(10,times = length(wspolczynniki))
cords <- list(x,y)
names(cords) <- c("x","y")
text(cords,labels = wspolczynniki,cex = 0.8)


# Zmutowana gameta --------------------------------------------------------
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
# skrypt do rysowania histogramu dla gamety
hist(replicate(100000,zmutowana.gameta(96,3.2)),freq = F,breaks = 15, main = "Liczba mutacji w gamecie przy liczbie mutacji 96 i wspó³czynniku 3.2",xlab = "Liczba mutacji",xlim = c(40,55),xaxp = c(40,55,15))

#skrypt porównuj¹cy iloœæ mutacji w gamecie w zale¿noœci od wsló³czynników oraz mutacji u rodzica
wspolczynniki <- c(1.2,3.2)
mut <- c(10,30,50,95)
kolory <- rainbow(length(wspolczynniki)-1)
kolory.plus <- c("black",kolory)
for(k in 1:length(mut)){
  for(j in 1:length(wspolczynniki)){
    replikaty <- replicate(100000,zmutowana.gameta(mut[k],wspolczynniki[j]))
    x <- sort(unique(replikaty))
    y <- c()
    for(i in 1:length(x)){
      y[i] <- sum(x[i] == replikaty)/length(replikaty)
    }
    if(j == 1 & k == 1){
      plot(x,y,type = "l",ylab = "Frekwencja",xlab = "Liczba mutacji w gamecie",xaxp = c(0,100,100),ylim = c(0,0.6),lwd = 2,xlim=c(0,100),main = "Funkcja - zmutowana gameta")
    } else if(k == 1){
      lines(x,y,type = "l",col = kolory[j-1],lwd=2)
    } else{
      lines(x,y,type = "l",col = kolory.plus[j],lwd=2,lty = k)
    }
  }
}
legend(52,0.3,title = "Wspó³czynniki",legend = wspolczynniki,col = kolory.plus,text.font = 4,pch = 15)
legend(73,0.3,title = "Liczba mutacji rodzica",legend = mut,text.font = 4,lty = 1:length(mut),lwd = 2)
