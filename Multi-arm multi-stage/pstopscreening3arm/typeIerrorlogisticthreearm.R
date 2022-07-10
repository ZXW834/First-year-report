load("nulllogisticvague.RData")

ntrials = 10000
ns = list(seq(25,225,25),seq(45,225,45),seq(75,225,75)) #group size of 15 for a maximum sample size of 150
null.response.probs = 0.12


pi.star=0.12
pess=2

model="ibb"
max.deviation=3
max.ar=0.8
rand.type=c("Coin","Urn")

PP.bound.grid=seq(0.99,0.999,0.0001)



library(parallel)

#cl <- makeCluster(40)

res={}
pstop={}
pstopscreening=function(ns,PP.bound.grid,res){
  Type.I.errorgridlist={}
  for (i in 1:3){
    Type.I.errorgrid={}
    for (z in 1:length(PP.bound.grid)){
      bounds= c(rep(PP.bound.grid[z],length(ns[[i]])-1),PP.bound.grid[z])  
      Type.I.Error = data.frame(rep(0,length(res)))
      Type.I.Error1=sapply(res[[i]],function(x) max(x[,1] > bounds | x[,1] < (1-bounds)) == 1)
      Type.I.Error2=sapply(res[[i]],function(x) max(x[,2] > bounds | x[,2] < (1-bounds)) == 1)
      typeIerrororder=c(max(mean(Type.I.Error1),mean(Type.I.Error2)),min(mean(Type.I.Error1),mean(Type.I.Error2)))
      Type.I.errorgrid=rbind(Type.I.errorgrid,typeIerrororder)
    }
    Type.I.errorgridlist=cbind(Type.I.errorgridlist,Type.I.errorgrid)
  }
  return(Type.I.errorgridlist)
}
pstop={}
Type.I.errorgrid=pstopscreening(ns=ns,PP.bound.grid = PP.bound.grid,res = res1)
Type.I.errorgrid2=pstopscreening(ns=ns,PP.bound.grid = PP.bound.grid,res = res2)

pstop25={}
pstop45={}
pstop75={}
for (i in 1:dim(Type.I.errorgrid)[1]){
  pstop25=c(pstop25,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,1]<=0.05)&&(Type.I.errorgrid[i,2]<=0.025))])
  pstop45=c(pstop45,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,3]<=0.05)&&(Type.I.errorgrid[i,4]<=0.025))])
  pstop75=c(pstop75,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,5]<=0.05)&&(Type.I.errorgrid[i,6]<=0.025))])
}

pstop=c(pstop,min(pstop25),min(pstop45),min(pstop75))

pstop25={}
pstop45={}
pstop75={}
for (i in 1:dim(Type.I.errorgrid)[1]){
  pstop25=c(pstop25,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,1]<=0.05)&&(Type.I.errorgrid2[i,2]<=0.025))])
  pstop45=c(pstop45,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,3]<=0.05)&&(Type.I.errorgrid2[i,4]<=0.025))])
  pstop75=c(pstop75,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,5]<=0.05)&&(Type.I.errorgrid2[i,6]<=0.025))])
}
pstop=c(pstop,min(pstop25),min(pstop45),min(pstop75))

names(pstop)<-c("25CoinV","45CoinV","75CoinV","25UrnV","45UrnV","75UrnV")
pstop


load("nulllogisticinfo.RData")

pstop={}
Type.I.errorgrid=pstopscreening(ns=ns,PP.bound.grid = PP.bound.grid,res = res1)
Type.I.errorgrid2=pstopscreening(ns=ns,PP.bound.grid = PP.bound.grid,res = res2)

pstop25={}
pstop45={}
pstop75={}
for (i in 1:dim(Type.I.errorgrid)[1]){
  pstop25=c(pstop25,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,1]<=0.05)&&(Type.I.errorgrid[i,2]<=0.025))])
  pstop45=c(pstop45,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,3]<=0.05)&&(Type.I.errorgrid[i,4]<=0.025))])
  pstop75=c(pstop75,PP.bound.grid[i*as.numeric((Type.I.errorgrid[i,5]<=0.05)&&(Type.I.errorgrid[i,6]<=0.025))])
}

pstop=c(pstop,min(pstop25),min(pstop45),min(pstop75))

pstop25={}
pstop45={}
pstop75={}
for (i in 1:dim(Type.I.errorgrid)[1]){
  pstop25=c(pstop25,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,1]<=0.05)&&(Type.I.errorgrid2[i,2]<=0.025))])
  pstop45=c(pstop45,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,3]<=0.05)&&(Type.I.errorgrid2[i,4]<=0.025))])
  pstop75=c(pstop75,PP.bound.grid[i*as.numeric((Type.I.errorgrid2[i,5]<=0.05)&&(Type.I.errorgrid2[i,6]<=0.025))])
}
pstop=c(pstop,min(pstop25),min(pstop45),min(pstop75))

names(pstop)<-c("25CoinL","45CoinL","75CoinL","25UrnL","45UrnL","75UrnL")
pstop