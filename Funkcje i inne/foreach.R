for(i in 1:4){
  g <- c(1:1000000)*i
  print(g)
}

library(foreach)
library(doParallel)
base <- 2
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
g <- c()
foreach(exponent = 2:4,.combine = c) %dopar%
  for(k in 5:10){
    k^exponent
  }
stopImplicitCluster()


install.packages("iterators", repos="http://R-Forge.R-project.org")
install.packages("foreach", repos="http://R-Forge.R-project.org")
install.packages("doParallel", repos="http://R-Forge.R-project.org")

# Przyk³ad 1 - ile czasu zajmuje krótkie for -----------------------------

iters<-10
ls<-vector('list',length=iters)
strt<-Sys.time()
for(i in 1:iters){
  cat(i,'\n')
  to.ls<-rnorm(1e6)
  to.ls<-summary(to.ls)
  ls[[i]]<-to.ls
}
print(Sys.time()-strt)

# Przyk³ad 2 - d³u¿sze for razem z kodem na fajny wykres-------------------

iters<-seq(10,100,by=10)
times<-numeric(length(iters))
for(val in 1:length(iters)){
  cat(val,' of ', length(iters),'\n')
  to.iter<-iters[val]
  ls<-vector('list',length=to.iter)
  strt<-Sys.time()
  for(i in 1:to.iter){
    cat(i,'\n')
    to.ls<-rnorm(1e6)
    to.ls<-summary(to.ls)
    ls[[i]]<-to.ls
  }
  times[val]<-Sys.time()-strt
}

library(ggplot2)

to.plo<-data.frame(iters,times)
ggplot(to.plo,aes(x=iters,y=times)) + 
  geom_point() + #dodaje punkty
  geom_smooth() + #linia dopasowania
  theme_bw() + #t³o wykresu
  scale_x_continuous('No. of loop iterations') + 
  scale_y_continuous ('Time in seconds')

mod<-lm(times~iters)
predict(mod,newdata=data.frame(iters=1e4))/60 #w minutach

# Przyk³ad 3 - proste foreach ---------------------------------------------
library(foreach)
library(doParallel)
iters <- 1e4
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
strt <- Sys.time()
writeLines(c(""),"log.txt")
ls <- foreach(i=1:iters,.combine = c) %dopar% {
  sink("log.txt",append = T)
  cat(paste("Iteration",i))
  to.ls <- rnorm(1e6)
  to.ls <- summary(to.ls)
  print(i)
  sink()
  to.ls
}
print(Sys.time()-strt)
stopCluster(cl)


# foreach z for w œrodku --------------------------------------------------
library(foreach)
library(doParallel)
iters <- c(1:100)
vec <- rep(1,times = 100)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
writeLines("","log.txt")
invisible(foreach(i=iters) %dopar% {
  sink("log.txt",append = T)
  to.ls <- c()
  to.ls[1] <- vec[1]
  for(j in 2:100){
    to.ls[j] <- to.ls[j-1] + vec[j]
  }
  to.ls <- paste(to.ls,collapse = " ")
  iteration <- paste("Iteration",i,sep = "")
  cat(paste(iteration,to.ls,"\n"))
  sink()
})
stopCluster(cl)

d <- read.table("log.txt",sep=" ",nrows = 1,skip = 1)
d2 <- as.numeric(d[1,2:(ncol(d)-1)])

# do instalowania "foreach" na Linuxie
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("foreach")
library(foreach)
