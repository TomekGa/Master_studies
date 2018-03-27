rozrzut <- read.table("choice_range1.1rep1mat3sconstanthconstant.txt", sep=" ", header = F)
mati <- matrix(c(1:50),nrow = 50,ncol = 1)
max.hist <- max(rozrzut)
bre <- seq(-0.5,18.5,by=1)

jpeg("choice_range1.1mat3sconstanthconstant.jpg",width = 500,height = 8000)

par(mfrow=c(50,1))

for(i in nrow(rozrzut):1){
  hist(x = as.numeric(rozrzut[i,]),col = "red",main = paste("Generation",i,sep = ""),xlim = c(0,18),xaxp=c(0,18,18),breaks = bre,ylim = c(0,500),
       xlab = "Number of possible homozygotes",ylab="Number of sperm cells")
}
dev.off()