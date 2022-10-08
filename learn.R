

#please first install packages "bnlearn" and "gRbase"

library(bnlearn)
library(gRbase)

##****************************************************************            
#==========            Directed Graphical Model    =============
#****************************************************************            

DAG<-function(simulate,burnin,Break,type1){

 data(learning.test)
 X<-learning.test[,c(3,6,1,2,4,5)]     # set ordering
 n<-nrow(X); K<-ncol(X)

 colnames(X)=1:K  
 Dag<-model2network("[1][2][3][4|3][5|1:3][6|2:4]")
 plot(Dag)
 B<-amat(Dag)
 score(Dag,X,type= type1)

#****************************************************************            
#      ==========    Simulation of the network    =========
#****************************************************************            

  E<-diag(rep(0,K))
  D.E<-c(); S.E<-c(); M.E<-c()

  R1<-0; sim<-0; 
  start_time <- Sys.time() 

  for(R in 1:simulate){

    print(R)

#****************************************************************            
#                 full conditional of E
#****************************************************************            

   for(i in 1:(K-1)){
     for(j in (i+1):K){
       E1<-E
       E1[i,j]<-0
       e<-empty.graph(as.character(1:K))
       amat(e)<-E1
       s0<-score(e,X,type= type1)    
       E1[i,j]<-1
       e<-empty.graph(as.character(1:K))
       amat(e)<-E1
       s1<-score(e,X,type= type1)
       ss<-abs(max(s0,s1))
       s0<-s0+ss; s1<-s1+ss; 
       p0<-exp(s0); p1<-exp(s1)
       E[i,j]<-sample(c(0,1),1,prob=c(p0,p1))
     }
   }      

#**************************************************************************
#                       Plot of convergence
#**************************************************************************

  S.E<-rbind(S.E,as.vector(E))
  M.E<-rbind(M.E,apply(S.E,2,mean))
 par(mfrow=c(3,5))
   plot(seq(1,R,1),M.E[,7], type="l",lwd=1,xlab="Iteration",ylab=expression(p[12]),main=expression(e[12]==0),cex.axis=.8)

   plot(seq(1,R,1),M.E[,13], type="l",lwd=1,xlab="Iteration",ylab=expression(p[13]),main=expression(e[13]==0),cex.axis=.8,col=1,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,14], type="l",lwd=1,xlab="Iteration",ylab=expression(p[23]),main=expression(e[23]==0),cex.axis=.8,col=2,ylim=c(0,1))

   plot(seq(1,R,1),M.E[,19], type="l",lwd=1,xlab="Iteration",ylab=expression(p[14]),main=expression(e[14]==0),cex.axis=.8,col=1,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,20], type="l",lwd=1,xlab="Iteration",ylab=expression(p[24]),main=expression(e[24]==0),cex.axis=.8,col=2,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,21], type="l",lwd=1,xlab="Iteration",ylab=expression(p[34]),main=expression(e[34]==1),cex.axis=.8,col=3,ylim=c(0,1))

   plot(seq(1,R,1),M.E[,25], type="l",lwd=1,xlab="Iteration",ylab=expression(p[15]),main=expression(e[15]==1),cex.axis=.8,col=1,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,26], type="l",lwd=1,xlab="Iteration",ylab=expression(p[25]),main=expression(e[25]==0),cex.axis=.8,col=2,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,27], type="l",lwd=1,xlab="Iteration",ylab=expression(p[35]),main=expression(e[35]==1),cex.axis=.8,col=3,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,28], type="l",lwd=1,xlab="Iteration",ylab=expression(p[45]),main=expression(e[45]==0),cex.axis=.8,col=4,ylim=c(0,1))

   plot(seq(1,R,1),M.E[,31], type="l",lwd=1,xlab="Iteration",ylab=expression(p[16]),main=expression(e[16]==0),cex.axis=.8,col=1,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,32], type="l",lwd=1,xlab="Iteration",ylab=expression(p[26]),main=expression(e[26]==1),cex.axis=.8,col=2,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,33], type="l",lwd=1,xlab="Iteration",ylab=expression(p[36]),main=expression(e[36]==0),cex.axis=.8,col=3,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,34], type="l",lwd=1,xlab="Iteration",ylab=expression(p[46]),main=expression(e[46]==1),cex.axis=.8,col=4,ylim=c(0,1))
   plot(seq(1,R,1),M.E[,35], type="l",lwd=1,xlab="Iteration",ylab=expression(p[56]),main=expression(e[56]==0),cex.axis=.8,col=5,ylim=c(0,1))

#****************************************************************            
#                             Save
#****************************************************************            

   if(R>burnin){
	  R1<-R1+1
      if(R1==Break){
	  sim<-sim+1		

       D.E<-rbind(D.E,as.vector(E))

#****************************************************************
       R1<-0

            }  #End if R1==break  
        }      #End if R>burnin
  }            #End for R.simulate

#****************************************************************
#                **************************
#                *     R E S U L T        *
#                **************************
#****************************************************************

print("E.P")
E.P<-matrix(apply(D.E,2,mean),K,K); print(round(E.P,2)); 

end_time <- Sys.time()
time<-end_time - start_time
cat("time====>",time,"\n")

#****************************************************************
print("*****************  R E S U L T  learn******************")
#****************************************************************

for(cut in seq(.4,.8,.05)){
  correct00<-0
  correct11<-0
  additional<-0
  missing<-0
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      if(B[i,j]==1 & E.P[i,j]>cut) correct11<-correct11+1
      if(B[i,j]==0 & E.P[i,j]<cut) correct00<-correct00+1
      if(B[i,j]==1 & E.P[i,j]<cut) missing<-missing+1
      if(B[i,j]==0 & E.P[i,j]>cut) additional<-additional+1
    }
  }
  totall<-K*(K-1)/2

  cat("cut point   ========================>>  ",cut,"\n")
  cat(" correct00====>>",correct00,"\n","correct11====>>",correct11,"\n",
      "percntage===>>",(correct00+correct11)/totall,"\n")
  cat(" PPV========>>",correct11/(correct11+additional),"\n",
     "TPR========>>",correct11/(correct11+missing),"\n")
  cat(" additional===>>",additional,"\n",
     "missing======>>",missing,"\n")
}

#********************************************* 
     
}   # END FUNCTION

#********************************************* 

DAG(simulate=200,burnin=100,Break=2,type1="k2")






