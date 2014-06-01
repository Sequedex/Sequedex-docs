# ======================================================================= 
#  Analysis protocol for species idenfitication in a metagenome
#  Nick Hengartner, June 6, 2013
#  Demo uses Ben's compiled dataset
#------------------------------------------------------------------------
#
#  step 1: clean the raw profiles by removing profiles that either 
#          have less than 1000 hits (mainly divergent organisms)
#          or profiles that have over 7.5% of reads that are not monophyletic
#          (possibly contaminated organisms)
#  step 2: cluster the remainin profiles into 800 groups
#  step 3: use the lasso to model the normalized metagenomic profile as a linear
#          combination of cluster profiles.  Use only nodes that are at depth 4 or less from
#          the leaves
#  step 4: extract the top 100 clusters (use RSS to assess if this number is too small)
#          that list is an upper bound of the possible profiles
#  step 5: Use traditional model selection approach to select significant profiles (use k=1.25
#          be more inclusive than AIC)
#-------------------------------------------------------------------------
#  Initialize functions, add packages

require(ape)
require(lars)

#-------------------------------------------------------------------------
#  Read in the reference profiles and metagenomic samples

rm(list=ls())
#setwd("~/LANL/Projects/Forensic/DemoSequedex")
source("sequestat4.r")
tree.file <- "Life2550.nexus"

#  input reference data and put them into the proper structure
syn  <- Read.tree(tree.file,colors=TRUE,sat=1)
tree <- syn
dat  <- read.table("syn.Life.who",header=F,sep="\t")  
lab  <- scan("syn.Life.lbl",what="",sep="\t")
lab  <- lab[-length(lab)]      # HACK  last entry is empty.  Remove.
colnames(dat) <- lab
nod.lab <- c(paste("n000",0:9,sep=""),
             paste("n00",10:99,sep=""),
             paste("n0",100:999,sep=""),
             paste("n",1000:(dim(dat)[1]-1),sep=""))
rownames(dat) <- nod.lab

syn$data <- dat
syn$node.label <- nod.lab

#  input metagenomic samples to be analyzed in profiles variable
dat=read.table("hmb.Life.who")
#  remove profiles with extreme low counts

metaG <- dat  
#idx <- apply(metaG,2,sum) > 1000
#metaG <- metaG[,idx]

#-------------------------------------------------------------------------
#  Step 1:  clean the reference profiles

all.paths <- All.paths(tree)
n.profiles <- dim(syn$data)[2]

sums.phylo <- Monophyl.sum(syn)
max.phylo <- apply(sums.phylo,2,which.max)
profil.phylo <- vector("list",n.profiles)
resid.phylo <- syn$data
for ( k in 1:n.profiles){
  profil.phylo[[k]] <- resid.phylo[all.paths[[max.phylo[k]]]$node,k]
  names(profil.phylo[[k]]) <- all.paths[[max.phylo[k]]]$node
  resid.phylo[all.paths[[max.phylo[k]]]$node,k] <- 0
}
names(profil.phylo) <- colnames(syn$data)

tree$data <- resid.phylo
err.phylo <- Monophyl.sum(tree)
err.ratio <- apply(resid.phylo,2,sum)/apply(syn$data,2,sum)
idx.bad   <- err.ratio > .075

tt <- table(apply(resid.phylo,2,sum))
tree <- syn
tree$data <- syn$data[,!idx.bad]


#-------------------------------------------------------------------------
#  step 2: cluster the remaining profiles into 800 groups

n.nodes  <- dim(tree$data)[1]
n.profil <- dim(tree$data)[2]
PP <- tree$data/matrix(apply(tree$data,2,sum),n.nodes,n.profil,byrow=T)
HH <- sqrt(PP)

set.seed(12345)
CC2 <- kmeans(t(HH),1000)
cc.name <- split(names(CC2$cluster),CC2$cluster)

#######################
#  build cluster archetype
n.cluster <- length(cc.name)
cP <- matrix(0,n.nodes,n.cluster)
for ( k in 1:n.cluster ){
  if ( CC2$size[k] == 1 ) cP[,k] <- tree$data[,CC2$cluster == k]
  if ( CC2$size[k] > 1 )  cP[,k] <- apply(tree$data[,CC2$cluster == k],1,sum)
}

#--------------------------------------------------------------------
#  step 3: use the lasso to model the normalized metagenomic profile as a linear

# tree$data <- cP
ndepth <- node.depth(tree)
ndepth <- ndepth[ ndepth > 1]

PP <- cP/matrix(apply(cP,2,sum),n.nodes,dim(cP)[2],byrow=T)
PP <- PP[ndepth < 251,]   # HACK  Focus only on outher shell --- discard inner nodes

#=========================================================
#  apply algorithm to a selected metagenome

meta.idx <- 1

Y0 <- metaG[,meta.idx]/sum(metaG[,meta.idx])
Y0 <- Y0[ndepth < 251 ]
fit0 <- lars(PP,Y0,type="lasso",intercept=F,normalize=F,max.steps=300,use.Gram=F)

#--------------------------------------------------------
#  step 4: extract the top 100 clusters

nvar <- 200
idx <- max((1:length(fit0$df))[fit0$df <= nvar])
rss <- fit0$RSS
bb <- fit0$beta[idx,]
cc <- bb[bb != 0]
big.list <- cc.name[bb != 0]
big.list <- big.list[sort.list(bb[bb != 0],decreasing=T)]

#-----------------------------------------------------------
#  step 5: Use traditional model selection approach to select significant profiles

vnames <- paste("V",1:dim(PP)[2],sep="")
colnames(PP) <- vnames
names(cc.name) <- vnames
mnames <- vnames[bb != 0]
mmodel0 <- formula(paste("Y0 ~ 0 + ",mnames[1]))
mmodel1 <- formula(paste("Y0 ~ ",paste(mnames,collapse=" + ")))
fit1 <- lm( mmodel0 ,data=as.data.frame(PP))
fit2 <- step(fit1, mmodel1, data=as.data.frame(PP),k=2)

sig.list <- cc.name[(names(fit2$coef))]
sig.list <- sig.list[sort.list(fit2$coef,decreasing=T)]
names(sig.list) <- paste(names(sig.list),sort(round(fit2$coef,5),decreasing=T))
