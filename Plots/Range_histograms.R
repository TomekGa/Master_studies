rozrzut <- read.table("choice_range1.1rep1mat1sconstanthconstant.txt", sep=" ", header = F)
mati <- matrix(c(1:50),nrow = 50,ncol = 1)
max.hist <- max(rozrzut)
bre <- seq(-0.5,23.5,by=1)

jpeg("choice_range1.1mat1sconstanthconstant.jpg",width = 500,height = 15000)

par(mfrow=c(50,1))

for(i in nrow(rozrzut):1){
  hist(x = as.numeric(rozrzut[i,]),col = "red",main = paste("Generation",i,sep = ""),xlim = c(0,23),xaxp=c(0,23,23),breaks = bre,ylim = c(0,200),
       xlab = "Number of mutations",ylab="Number of sperm cells",cex.lab=1.5)
}
dev.off()