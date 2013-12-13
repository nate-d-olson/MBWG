# determination of full length 16S individual gene copy sequences
# using maximum likelihood
# written by: Nate Olson, statistical analysis method developed by Steve Lund

library(ggplot2)
library(reshape2)
Dat<-read.csv("~/Documents/CCQM/ecoli_full_length_analysis/ecoli_copy_proportions.csv")
dat_df <- as.data.frame(Dat[-1])
colnames(dat_df) <- c("Variant String","LGC_454.1","LGC_454.2","LGC_454.3","NMIA_454","NIST_Sanger","LGC_Sanger")
rownames(Dat)<-Dat[,2]
Dat<-Dat[,c(-1,-2)]
dat<-rowSums(as.matrix(Dat))

pats<-cbind(rep(0:1,each=2),c(0,1,0,1))

P.vec<-function(p.switch,allele.ind){
  p.vec<-NULL
  for(i in 1:nrow(pats)){
    p.vec<-c(p.vec,(1-p.switch)*mean(allele.ind==i)+
               p.switch*mean(pats[allele.ind,2]==pats[i,2]))
  }
  p.vec
}

like<-function(p.switch,allele.ind,dat) 
  -dmultinom(dat,prob=P.vec(p.switch,allele.ind),log=TRUE)


like(.005,c(1,1,1,2,2,2,4),dat)

res<-NULL
for(i in 1:4){
  print(i)
  for( j in 1:i ){
    for( k in 1:j ){
      for( l in 1:k ){
        for( m in 1:l ){
          for( n in 1:m ){
            for( o in 1:n){
              t.ind<-c(i,j,k,l,m,n,o)
              fit<-optimize(like,c(0,1),allele.ind=t.ind,dat=dat)
              res<-rbind(res,c(t.ind,fit$minimum,fit$objective))  
            } } } } } }}

dat_df$Predicted <- c(3,3,0,1)
dat_dfm <- melt(dat_df)
dat_dfc <- dcast(dat_dfm, variable~., sum)
colnames(dat_dfc) <- c("variable", "sample_size")
ggplot(dat_dfm) + geom_bar(aes(x = variable, y= value , fill= string), stat = "identity", position = "fill") +
  geom_text(data = dat_dfc, aes(x = variable, y = 1.05, label = sample_size)) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Dataset", y = "Proportion", fill = "Variant String")


res_df <- as.data.frame(res)
colnames(res_df) <- c(paste(rep("copy_", 7),1:7, sep = ""), "chimera","likelihood")
optimum <- min(res_df$likelihood)
ggplot(res_df) + geom_point(aes(x = likelihood, y = chimera)) + 
  geom_point(data = res_df[res_df$likelihood == optimum,], aes(x = likelihood, y = chimera), color = "red", size = 4)  + 
  scale_x_log10() + scale_y_log10() +
  labs(x = "Log Likelihood", y = "Log Chimera Probability")
