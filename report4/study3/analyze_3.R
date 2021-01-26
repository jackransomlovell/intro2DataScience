
# load packages -----------------------------------------------------------
if(!require(psych)) install.packages("psych"); require(psych)
if(!require(pracma)) install.packages("pracma"); require(pracma)
if(!require(lattice)) install.packages("lattice"); require(lattice)
if(!require(grid)) install.packages("grid"); require(grid)
if(!require(expm)) install.packages("expm"); require(expm)
if(!require(ppcor)) install.packages("ppcor"); require(ppcor)

# load data ---------------------------------------------------------------
# rating data
tlong <- read.csv("transitions_long.csv")
flong <- read.csv("frequencies_long.csv")
rwide <- read.csv("ratings_wide.csv")

# experience data
freqs <- read.csv("ground_frequencies.csv")
gfreqs <- freqs[,2]
names(gfreqs)<-tolower(freqs[,1])

gtrans <- read.csv("ground_transition_odds.csv")
gtrans <- gtrans[,2:dim(gtrans)[2]]
gtrans <- as.matrix(gtrans)
rownames(gtrans)<-tolower(colnames(gtrans))
colnames(gtrans)<-tolower(colnames(gtrans))

gcooc <- read.csv("ground_co-occurence_odds.csv")
gcooc <- gcooc[,2:dim(gcooc)[2]]
gcooc <- as.matrix(gcooc)
rownames(gcooc)<-tolower(colnames(gcooc))
colnames(gcooc)<-tolower(colnames(gcooc))

# reliability -------------------------------------------------------------

# frequency
freqs <- unname(t(rwide[,329:346]))
cormat <- cor(freqs)
k <- dim(cormat)[1]
diag(cormat)<-0
mr <- mean(squareform(cormat)) # average interrater
std.alpha <- k*mr/(1+(k-1)*mr) # interrater alpha

# frequency (attention passed)
freqs <- unname(t(rwide[rwide$attnpassed==1,329:346]))
cormat <- cor(freqs)
k <- dim(cormat)[1]
diag(cormat)<-0
mr <- mean(squareform(cormat)) # average interrater
std.alpha <- k*mr/(1+(k-1)*mr) # interrater alpha

# transitions
trans <- unname(t(rwide[,5:328]))
cormat <- cor(trans)
k <- dim(cormat)[1]
diag(cormat)<-0
mr <- mean(squareform(cormat)) # average interrater
std.alpha <- k*mr/(1+(k-1)*mr) # interrater alpha

# transitions (attention passed)
trans <- unname(t(rwide[rwide$attnpassed==1,5:328]))
cormat <- cor(trans)
k <- dim(cormat)[1]
diag(cormat)<-0
mr <- mean(squareform(cormat)) # average interrater
std.alpha <- k*mr/(1+(k-1)*mr) # interrater alpha

# test significance of transition inter-rater r
trans <- unname(t(rwide[,5:328]))
tcor <- mean(1-squareform(1-cor(trans)))
# (approach 1 - permuting ratings)
tperm <- replicate(1000,mean(1-squareform(1-cor(apply(trans,2,function(x) sample(x))))))
mean(tcor<abs(tperm))
# (approach 2 - permuting rows and columns)
tkey <- tlong[1:324,5:6]
indmark<-tapply(tlong$response,tlong$id,function(x) tapply(x,tkey,mean))
permcol <- function(mat){
  n <- dim(mat)[1]
  sel <- sample(n)
  return(mat[sel,sel])
}
set.seed(1)
tperm <- replicate(1000, mean(1-squareform(1-cor(do.call(cbind,lapply(indmark,function(x) as.vector(permcol(x))))))))
mean(tcor<abs(tperm))

# frequency: intuition vs. ground truth -----------------------------------

ifreqs <- 100-colMeans(rwide[,329:346])
names(ifreqs)<-tolower(as.character(flong$state)[1:length(ifreqs)])
ifreqs <- ifreqs[names(gfreqs)]
cor(ifreqs,gfreqs,method="spearman") # correlation

tiff("ifreq_gfreq_scatter.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3,3,1,1),mgp =c(2,.5,0))
plot(ifreqs,gfreqs*100,xlab="Rated frequencies (% of time)",ylab="Experienced frequencies (% of time)",
     pch=20)
yshift <- rep(-.75,length(ifreqs))
xshift <- rep(0,length(ifreqs))
xshift[9] <- -1
xshift[15] <- 1
xshift[13] <- -1
xshift[12] <- .5
yshift[15] <- 1
yshift[13] <- 1
yshift[16] <- -1.2
text(ifreqs+xshift,gfreqs*100+yshift,names(ifreqs))
text(39,34,"rating = experience")
text(60,17,"ρ = .75")
abline(1:100,1:100)
abline(lm(gfreqs*100~ifreqs),lty=2)
dev.off()


# individual raters
ifreqs <- 100-rwide[,329:346]
names(ifreqs)<-tolower(as.character(flong$state)[1:length(ifreqs)])
ifreqs <- ifreqs[names(gfreqs)]
irs <- cor(t(ifreqs),gfreqs,method="spearman")
set.seed(1)
boots <- quantile(replicate(10000,mean(sample(irs,size = length(irs),T))),c(.025,.975))
tiff("indivifreq_gfreq.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3,3,1,1),mgp =c(2,.5,0))
hist(irs,10,xlab="Individual correlations (ρ) with average experience",main=NULL,
     col="grey",border="darkgrey",ylim=c(0,25),xlim=c(-1,1))
abline(v=mean(irs))
abline(v=boots,lty=2)
dev.off()

# transitions: ratings and experienced ------------------------------------

color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}
steps <- c("blue4", "cyan", "white", "yellow", "red4")
pal <- color.palette(steps, c(50,1,1,50), space="rgb")

# rating matrix
tlong$state1<-tolower(tlong$state1)
tlong$state2<-tolower(tlong$state2)
intrans <- tapply(tlong$response,tlong[,5:6],mean)
intrans <- intrans[names(sort(gfreqs)),names(sort(gfreqs))]
tiff("transition_ratings.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
levelplot(t(intrans),main=list(label="Intuitive transition likelihoods",cex=1.5),
          ylab=list(label="From",cex=1.5),xlab=list(label="To",cex=1.5),at=seq(0,100,1), 
          scales=list(x=list(rot=45,cex=1.2),y=list(cex=1.2)),col.regions=pal,colorkey=list(labels=list(cex=1.2)))
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text("%",.4,-.05,1,1,gp=gpar(fontsize=18))
trellis.unfocus()
dev.off()


# experienced matrix
colnames(gtrans)<-gsub("\\.","-",colnames(gtrans))
rownames(gtrans)<-gsub("\\.","-",rownames(gtrans))
gtrans <- gtrans[names(sort(gfreqs)),names(sort(gfreqs))]
colnames(gcooc)<-gsub("\\.","-",colnames(gcooc))
rownames(gcooc)<-gsub("\\.","-",rownames(gcooc))
gcooc <- gcooc[names(sort(gfreqs)),names(sort(gfreqs))]
tiff("transition_experience.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
levelplot(t(log(gtrans)),main=list(label="Experienced transition log odds",cex=1.5),
          ylab=list(label="From",cex=1.5),xlab=list(label="To",cex=1.5),at=seq(-2.7,2.7,length.out = 100), 
          scales=list(x=list(rot=45,cex=1.2),y=list(cex=1.2)),col.regions=pal,colorkey=list(labels=list(cex=1.2)))
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
grid.text("log odds",1,-.05,1,1,gp=gpar(fontsize=18))
trellis.unfocus()
dev.off()

# rating-experience correlation
invec <- as.vector(intrans)
gvec <- as.vector(log(gtrans))
gcvec <- as.vector(log(gcooc))
cor(invec,gvec,method="spearman")
tiff("transition_scatter.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3.6,4,1,1),mgp =c(2.4,.8,0))
plot(invec,gvec,xlab="Mental model transition probability (%)",ylab="Experienced transition (log odds)",
     pch=20,xlim=c(10,100),cex.axis=1.6,cex.lab=2)
text(97,.75,"ρ = .79",cex=1.6)
abline(lm(gvec~invec),lty=2,lwd=3)
dev.off()


# individual rating-experience correlations
tkey <- tlong[1:324,5:6]
norder <- names(sort(gfreqs))
indmark<-tapply(tlong$response,tlong$id,function(x) tapply(x,tkey,mean)[norder,norder])
indcors <- unlist(lapply(indmark,function(x) cor(as.vector(x),gvec,method="spearman")))
co.indcors<- unlist(lapply(indmark,function(x) pcor(cbind(as.vector(x),gcvec,gvec),method="spearman")$estimate[1,3]))
mean(indcors)/sd(indcors)
set.seed(1)
boots <- quantile(replicate(10000,mean(sample(indcors,size = length(indcors),T))),c(.025,.975))
perms <- replicate(10000,mean(unlist(lapply(indmark,function(x) cor(as.vector(x),gvec[sample(length(gvec))],method="spearman")))))
mean(abs(perms)>mean(indcors))
tiff("indivitrans.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3.6,4,1,1),mgp =c(2.4,.8,0))
hist(indcors,10,xlab="Correlations (ρ) with experienced transitions",main=NULL,
     col="grey",border="grey",ylim=c(0,40),xlim=c(-.5,1),cex.axis=1.6,cex.lab=2)
abline(v=mean(indcors),lwd=3)
abline(v=boots,lty=2,lwd=3)
dev.off()

# cross-validation --------------------------------------------------------

# with respect to emotions
set.seed(1)
nfold <- 5
cvinds <- sample(rep(1:nfold,ceil(dim(intrans)[1]/nfold)),dim(intrans)[1])
rmse <- rep(NA,nfold)
rrmse <- rep(NA,nfold)
for (i in 1:nfold){
  sel <- cvinds == i
  y <- as.vector(intrans[!sel,!sel])
  x <- as.vector(log(gtrans)[!sel,!sel])
  cvfit <- lm(y~x)
  sintrans <- intrans[sel,sel]
  y<-as.vector(sintrans)
  rvec <- sample(sum(sel),sum(sel))
  yr<-as.vector(sintrans[rvec,rvec])
  x <- as.vector(log(gtrans)[sel,sel])
  rmse[i] <- sqrt(mean((y-predict(cvfit,data.frame(x=x)))^2))
  rrmse[i] <- sqrt(mean((yr-predict(cvfit,data.frame(x=x)))^2))
}
mean(rmse)
mean(rrmse)
mean(rmse)/diff(range(intrans))
mean(rrmse)/diff(range(intrans))


# with respect to participants
set.seed(1)
nfold <- 5
cvinds <- sample(rep(1:nfold,ceil(dim(rwide)[1]/nfold)),dim(rwide)[1])
mrmse <- rep(NA,nfold)
mrrmse <- rep(NA,nfold)
for (i in 1:nfold){
  sel <- cvinds == i
  yvars<-do.call(cbind,lapply(indmark[!sel],function(x) as.vector(x)))
  cvfit <- rowMeans(apply(yvars,2,function(y) lm(y~gvec)$coefficients))
  yvars <- do.call(cbind,lapply(indmark[sel],function(x) as.vector(x)))
  rvec <- sample(dim(intrans)[1],dim(intrans)[1])
  ryvars <- do.call(cbind,lapply(indmark[sel],function(x) as.vector(x[rvec,rvec])))
  mrmse[i] <- mean(apply(yvars,2,function(x) rmserr(x,gvec*cvfit[2]+cvfit[1])$rmse))
  mrrmse[i] <- mean(apply(ryvars,2,function(x) rmserr(x,gvec*cvfit[2]+cvfit[1])$rmse))
}
mean(mrmse)
mean(mrrmse)
mean(rmse)/diff(range(tlong$response))
mean(rrmse)/diff(range(tlong$response))


# egocentrism -------------------------------------------------------------

freqs <- unname(t(rwide[,329:346]))
cf <- cor(freqs)
cf <- fisherz(cf)
diag(cf)<-0
trans <- unname(t(rwide[,5:328]))
ct <- cor(trans)
ct <- fisherz(ct)
diag(ct)<-0
cor(squareform(cf),squareform(ct))

cperm <- function(cmat){
  n <- dim(cmat)[1]
  sel <- sample(n,n)
  cmat <- cmat[sel,sel]
  return(cmat)
}
set.seed(1)
permvals <- replicate(10000,cor(squareform(cperm(cf)),squareform(ct)))
mean(cor(squareform(cf),squareform(ct))<abs(permvals))



# markov stationary distribution test -------------------------------------

# item analysis
intransnorm <- intrans/rowSums(intrans)*100
stationary <- ((intransnorm/100)%^%(10))[1,]
cor(stationary,gfreqs,method = "spearman")
tiff("stationary_gfreq_scatter.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3,3,1,1),mgp =c(2,.5,0))
plot(stationary*100,gfreqs*100,xlab="Markov chain stationary distribution (% of time)",ylab="Experienced frequencies (% of time)",
     pch=20)
yshift <- rep(-.75,length(ifreqs))
xshift <- rep(0,length(ifreqs))
xshift[12]<-.02
xshift[7]<-.02
xshift[9]<--.02
text(stationary*100+xshift,gfreqs*100+yshift,names(ifreqs))
text(6,21,"ρ = .64")
abline(lm(gfreqs*100~I(stationary*100)),lty=2)
dev.off()

# individual participants
frats <- t(100-rwide[,329:346]) 
rownames(frats)<-tolower(as.character(flong$state)[1:18])
frats <- frats[names(gfreqs),]
indcors <-do.call(rbind,lapply(indmark,function(x) cor(((x/rowSums(x))%^%(10))[1,],frats,method = "spearman")))
dindcors <- diag(indcors) # to their own frequency ratings
# to the experienced frequencies
indacc <- unlist(lapply(1:length(indmark),function(x) cor(((indmark[[x]]/rowSums(indmark[[x]]))%^%(10))[1,],gfreqs,method = "spearman")))
# partial correlation
indpar <- unlist(lapply(1:length(indmark),function(x) pcor.test(((indmark[[x]]/rowSums(indmark[[x]]))%^%(10))[1,],gfreqs,frats[,x],method = "spearman")$estimate))

set.seed(1)
boots1 <- quantile(replicate(10000,mean(sample(dindcors,size = length(dindcors),T))),c(.025,.975))
boots2 <- quantile(replicate(10000,mean(sample(indacc,size = length(indacc),T))),c(.025,.975))
boots3 <- quantile(replicate(10000,mean(sample(indpar,size = length(indpar),T))),c(.025,.975))
boots4 <- quantile(replicate(10000,mean(sample(fisherz(indacc)-fisherz(indpar),size = length(indacc-indpar),T))),c(.025,.975))


tiff("indivistationary.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(3,3,1,1),mgp =c(2,.5,0))
hist(indacc,10,xlab="Individual stationary distribution correlations (ρ) with average experienced transitions",main=NULL,
     col="grey",border="darkgrey",ylim=c(0,15),xlim=c(-1,1))
abline(v=mean(indacc))
abline(v=boots2,lty=2)
dev.off()


# markov chain foresight comparison ---------------------------------------

markovchoose <- function(tpm,current){
  nstate <- dim(tpm)[1]
  currow <- tpm[current,]
  newstate <- sample(1:nstate,1,F,currow)
  return(newstate)
}


set.seed(1)

# simulate random experience walks
nwalk <- 10000
nstep <- 4
wstart <- sample(1:18,nwalk,replace = T)
tpm <- gtrans/rowSums(gtrans)
walks <- matrix(NA,nwalk,nstep)
for (i in 1:nwalk){
  curstate <- wstart[i]
  for (j in 1:nstep){
    curstate<-markovchoose(tpm,curstate)
    walks[i,j]<-curstate
  }
}

nsub <- dim(rwide)[1]
acc <- matrix(NA,nsub,nstep)
for (i in 1:nsub){
  tpm <- indmark[[i]]
  tpm <- tpm/rowSums(tpm)
  curacc <- matrix(NA,nwalk,nstep)
  for (j in 1:nwalk){
    curwalk <- rep(NA,nstep)
    curstate <- wstart[j]
    for (k in 1:nstep){
      curstate<-markovchoose(tpm,curstate)
      curwalk[k]<-curstate
    }
    curacc[j,] <- as.numeric(curwalk==walks[j,])
  }
  acc[i,]<-colMeans(curacc)
}
m <- colMeans(acc)
colMeans(acc>(1/18))


set.seed(1)
bootsteps <- apply(replicate(10000,colMeans(acc[sample(1:nsub,nsub,T),])),1,function(x) quantile(x,c(.025,.975)))


tiff("foresight.tif",width=20,height=20,units="cm",res=600,compression = "lzw")
par(mar=c(4,4,1,1),mgp =c(2.6,1,0))
plot(1:4,m,type="o",pch=20,lty=3,cex=3,lwd=3,
     xlim=c(0,4),ylim=c(.054,.068),cex.axis=1.6,cex.lab=2,
     xlab="Steps into future",ylab="Simulated accuracy proportion")
abline(h=1/18,lty=2,lwd=2)
abline(v=0,lty=2,lwd=2)
segments(1:4,bootsteps[1,],1:4,bootsteps[2,],lwd=3)
text(.1,.059,"Observed emotion",srt=90,cex=1.6)
text(.5,.0552,"Chance",cex=1.6)
dev.off()


# Valence analysis --------------------------------------------------------

valmod <- matrix(0,18,18)
valmod[1:9,1:9]<-1
valmod[10:18,10:18]<-1
  
valcor <- unlist(lapply(indmark,function(x) cor(as.vector(valmod),as.vector(x),method = "spearman")))
mean(valcor)/sd(valcor)
