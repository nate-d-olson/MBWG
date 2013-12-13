# determination of full length 16S individual gene copy sequences
# using maximum likelihood
# written by: Nate Olson, statistical analysis method developed by Steve Lund

Dat<-read.csv("copy_proportions.csv")
rownames(Dat)<-Dat[,1]
Dat<-Dat[,-1]
dat<-rowSums(as.matrix(Dat))
 
pats<-cbind(rep(0:1,each=4),c(0,0,1,1,0,0,1,1),rep(0:1,4))
 
P.vec<-function(p.switch,allele.ind){
                p.switch.short<-1-(1-p.switch)^13
                p.switch.long<-1-(1-p.switch)^131
                p.vec<-NULL
                for(i in 1:nrow(pats)){
                                p.vec<-c(p.vec,(1-p.switch.short)*(1-p.switch.long)*mean(allele.ind==i)+
                                (1-p.switch.short)*p.switch.long*mean(paste(pats[allele.ind,1],pats[allele.ind,2])==paste(pats[i,1],pats[i,2]))*mean(pats[allele.ind,3]==pats[i,3])+
                                p.switch.short*(1-p.switch.long)*mean(paste(pats[allele.ind,2],pats[allele.ind,3])==paste(pats[i,2],pats[i,3]))*mean(pats[allele.ind,1]==pats[i,1])+
                                p.switch.short*p.switch.long*mean(pats[allele.ind,1]==pats[i,1])*mean(pats[allele.ind,2]==pats[i,2])*mean(pats[allele.ind,3]==pats[i,3]))
                }
                p.vec
}
 
like<-function(p.switch,allele.ind,dat) -dmultinom(dat,prob=P.vec(p.switch,allele.ind),log=TRUE)
 
 
like(.005,c(2,2,3,4,4,6),dat)
 
res<-NULL
for(i in 1:8){
 print(i)
for( j in 1:i ){
  for( k in 1:j ){
   for( l in 1:k ){
    for( m in 1:l ){
     for( n in 1:m ){
                t.ind<-c(i,j,k,l,m,n)
                fit<-optimize(like,c(0,1),allele.ind=t.ind,dat=dat)
                res<-rbind(res,c(t.ind,fit$minimum,fit$objective))
} } } } } }
