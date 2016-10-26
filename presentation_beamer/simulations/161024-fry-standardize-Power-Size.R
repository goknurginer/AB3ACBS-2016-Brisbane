# 24 October 2016
library(limma)
library(ggplot2)

ngenes <- 100000
y <- matrix(rnorm(ngenes*6),ngenes,6)

design <- cbind(1,c(0,0,0,1,1,1))

s0 <- 0.3
d0 <- 4
sigma2 <- s0^2 * d0 / rchisq(ngenes, df=d0)
y <- y*sqrt(sigma2)

nsets <- ngenes/5 # 20000 sets
Truth <- rep(0,nsets)
Truth[sample(nsets,2000)] <- 1 # randomly select the expressed sets which is only the 10% of the sets
table(Truth) # 2000 of the sets are expressed
index <- list() # generate set indices. Each set has only 5 genes
i <- 1:5 # Add 2 to the second factor(samples 4:6) of the genes which falls into expressed sets
for(s in 1:nsets) {
	index[[s]] <- i
	if(Truth[s]) y[i,4:6] <- y[i,4:6] + 2
	i <- i+5
}

f.p <- fry(y,index=index,design=design,sort="n") # standardize = "posterior.sd" use default
ro <- mroast(y,index=index,design=design,nrot=999,sort="n")

# Check if fry and roast hold tehir sizes
p <- cbind(f=f.p$PValue,ro=ro$PValue)
colMeans(p<0.05)
colMeans(p<0.01)
colMeans(p<0.001)
colMeans(p<0.0001)

o <- order(f.p$PValue)
NFD.p <- cumsum(1-Truth[o])
o <- order(ro$PValue)
NFD.ro <- cumsum(1-Truth[o])

n <- 2000 # number of expressed sets
m <- max(NFD.ro[n],NFD.p[n])
plot(1:n,seq(0,m,len=n),type="n",xlab="Number selected",ylab="Number of false discoveries")
lines(1:n,NFD.p[1:n],col="red")
lines(1:n,NFD.ro[1:n],col="black",lty=2)
legend("topleft",lty=c(1,2),col=c("red","black"),legend=c("fry","roast"))


# ggplot
df.p <- data.frame(Method=rep("Fry", nsets), NFD=NFD.p[1:n])
df.ro <- data.frame(Method=rep("Roast", nsets), NFD=NFD.ro[1:n])
df.p <- df.p[order(df.p$NFD),]
df.ro <- df.ro[order(df.ro$NFD),]
df <- rbind(df.p,df.ro)

p <- qplot(c(1:20000,1:20000), NFD, data = df, geom="step", colour = Method)

setEPS()
postscript(file = "..//figures//power1.eps", width = 10, height = 8)
p + xlab("Number of Selected Sets") + 
  ylab("Number of False Discoveries") + 
  theme_minimal()
dev.off()

















