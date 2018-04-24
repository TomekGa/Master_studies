#Wykresy porownujace macierze odleglosci 

l.pokolen <- 50
inbreeding.loci.nr <- 7
for(i in 1:l.pokolen){
  odleglosci.name <- paste("odleglosci",i,".txt",sep = "")
  odleglosci.name2 <- paste("odleglosci",i,sep = "")
  assign(odleglosci.name2,read.table(odleglosci.name,sep = "\t",header = T,colClasses = "character"))
}

jpeg("mat5_macierze_odleglosci.jpg",width = 500,height = 15000)
par(mfrow=c(50,1))

for(i in 50:1){
  odleglosci.name2 <- paste("odleglosci",i,sep = "")
  macierz <- eval(parse(text = odleglosci.name2))
  for(j in 1:nrow(macierz)){
    wiersz <- macierz[j,]
    x <- c(0:inbreeding.loci.nr)
    y <- c()
    for(k in 1:length(x)){
      y[k] <- sum(wiersz == x[k])
    }
    if(j == 1){
      plot(x=x,y = y)
    } else {
      lines(x = x,y = y)
    }
  }
}
dev.off()