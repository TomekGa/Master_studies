#par.pok<-read.table(file="parametry-pokolenia.txt", sep="\t", header = T)
par1<-read.table(file="parametry1variable.txt", sep="\t", header = T)
par2<-read.table(file="parametry2variable.txt", sep="\t", header = T)
par3<-read.table(file="parametry3variable.txt", sep="\t", header = T)

par11<-read.table(file="parametry1constant.txt", sep="\t", header = T)
par22<-read.table(file="parametry2constant.txt", sep="\t", header = T)
par33<-read.table(file="parametry3constant.txt", sep="\t", header = T)

pokolenia<-c(0:(nrow(par1)-1))
nazwy_parametrow<-colnames(par1)

#viability
par(mfrow=c(2,2))
kolory <- rainbow(2)
#kolory.plus <- c("black",kolory)
for (p in 1:4){
  if(p<4){
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab=1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(0.3,2.5))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  } else {
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab = 1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(0.01,0.35))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  }
}


#next 4
par(mfrow=c(2,2))
kolory <- rainbow(2)
#kolory.plus <- c("black",kolory)
for (p in 5:8){
  if(p<7){
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab=1.5,lwd = 4,ylab = nazwy_parametrow[p])
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  } else if(p==7) {
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab = 1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(1.5,120))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  } else if(p==8) {
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab = 1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(0,60))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  }
}

#next 3
par(mfrow=c(2,2))
kolory <- rainbow(2)
#kolory.plus <- c("black",kolory)
for (p in 9:11){
  if(p==9){
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab=1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(0,225))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4) 
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  } else if(p==10) {
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab = 1.5,lwd = 4,ylab = nazwy_parametrow[p])
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4)
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
    lines(x=pokolenia, y=par11[,p], type="l",col = "black",lwd = 4,lty = 2)
    lines(x=pokolenia,y=par22[,p],type = "l",col = kolory[1],lwd=4,lty = 2)
    lines(x=pokolenia,y=par33[,p],type = "l",col = kolory[2],lwd=4,lty = 2)
  } else if(p==11) {
    plot(x=pokolenia, y=par1[,p], type="l", xlab="pokolenia",cex.lab = 1.5,lwd = 4,ylab = nazwy_parametrow[p],ylim = c(14,21))
    lines(x=pokolenia,y=par2[,p],type = "l",col = kolory[1],lwd=4)
    lines(x=pokolenia,y=par3[,p],type = "l",col = kolory[2],lwd=4)
  }
}
# par(op) # Leave the last plot
# op <- par(usr=c(0,1,0,1), # Reset the coordinates
#           xpd=NA)
# legend(-0.7,1.75,title = "mating system",legend = c("Males with the smallest nr of mutations","Random mating","Smallest nr of common mutations"),col = c("black",kolory),text.font = 4,pch = 15)

#b) Je?li chcemy narysowa? pojedynczo jaki? jeden wykres:
parametr<-"mean_viability" #tu nazwa parametru, dla kt?rego chcemy wykres
zakres<-c()
par(mfrow=c(1,1))
plot(x=pokolenia, y=par.[,parametr], type="l", xlab="pokolenia", ylab=parametr, ylim=zakres)8
