rep <- 1
file.name3 <- paste("replicate",rep,"plot.txt",sep="")
l.pokolen <- 3
Nmax <- 10
coordinatesTable<- read.table("plotrep1mat4sconstanthconstant.txt",sep = "\t",header = T,colClasses = "character")
x.plot.homo.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[1,])))
y.plot.homo.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[2,])))
x.plot.hetero.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[3,])))
y.plot.hetero.1 <- as.numeric(na.omit(as.numeric(coordinatesTable[4,])))
x.plot.homo.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[5,])))
y.plot.homo.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[6,])))
x.plot.hetero.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[7,])))
y.plot.hetero.0 <- as.numeric(na.omit(as.numeric(coordinatesTable[8,])))

x.plot.homo.1 <- x.plot.homo.1[-1]
y.plot.homo.1 <- y.plot.homo.1[-1]
x.plot.hetero.1 <- x.plot.hetero.1[-1]
y.plot.hetero.1 <- y.plot.hetero.1[-1]
x.plot.homo.0 <- x.plot.homo.0[-1]
y.plot.homo.0 <- y.plot.homo.0[-1]
x.plot.hetero.0 <- x.plot.hetero.0[-1]
y.plot.hetero.0 <- y.plot.hetero.0[-1]

yax <- ((1+l.pokolen)*Nmax)
opisy <- rep(1:10,times=l.pokolen+1)
pokolenia <- c(0:l.pokolen)
jpeg("rplot.jpg",width = 15000,height = 6000)
mar.default <- c(5,4,4,2) + 0.1
par(xpd = NA,oma=rep(0,times=4),mar=mar.default + c(4, 4, 4, 4))
plot(x = x.plot.hetero.1,y = y.plot.hetero.1,xaxp = c(0,6000,6000),yaxp=c(0,yax,yax),ylim = c(0,yax),pch=17,col="black",yaxt="n",xlab = "loci",ylab = "Number of individuals and generation (in red)",cex.lab = 3,cex.axis=3,cex=3,xaxs="i",yaxs="i")
axis(2, at=seq(0, yax, by=1), labels = FALSE)
text(par("usr")[1]-0.5,seq(1, yax, by=1),labels=opisy,pos=2,cex = 3)
text(par("usr")[1]-20,seq((Nmax/2), yax-(Nmax/2), by=Nmax),labels=pokolenia,pos=2,col="red",cex = 3)
points(x=x.plot.homo.1,y=y.plot.homo.1,pch=17,col="red",cex=3)
points(x=x.plot.hetero.0,y=y.plot.hetero.0,pch=19,col="black",cex=3)
points(x=x.plot.homo.0,y=y.plot.homo.0,pch=19,col="red",cex=3)
lines(x=c(3000.5,3000.5),y=c(0,yax))
for(i in 1:l.pokolen){
  lines(x=c(0,6000),y=c((i*Nmax)+0.5,y=c((i*Nmax)+0.5)))
}

dev.off()
