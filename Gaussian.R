

#please first install packages "bnlearn" and "gRbase"

library(bnlearn)
library(gRbase)

##****************************************************************            
#==========          Directed Graphical Model    =============
#****************************************************************            

DAG<-function(simulate,burnin,Break){

 data(gaussian.test)
 X<-gaussian.test
 Dag0= model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
 order<-topo_sort(amat(Dag0)) # gRbase package
 X<-X[order]                  # set ordering
 B<-amat(Dag0)[order,order]   # set ordering
 n<-nrow(X); K<-ncol(X)
 colnames(X)=1:K 
 colnames(B)=1:K 
 rownames(B)=1:K 
 e = empty.graph(as.character(1:K))
 amat(e)=B
 plot(e)

#****************************************************************            
#      ==========    Simulation of the network    =========
#****************************************************************            

  E<-diag(rep(0,K))
  D.E<-c(); S.E<-c(); M.E<-list()

  R1<-0; sim<-0; 
  start_time <- Sys.time() 

  for(R in 1:simulate){

    cat("simulation=",R,"\n")

#****************************************************************            
#                 full conditional of E
#****************************************************************            

   for(i in 1:(K-1)){
     for(j in (i+1):K){
       E1<-E
       E1[i,j]<-0
       e<-empty.graph(as.character(1:K))
       amat(e)<-E1
       s0<-score(e,X,type="bge")    
       E1[i,j]<-1
       e<-empty.graph(as.character(1:K))
       amat(e)<-E1
       s1<-score(e,X,type="bge")
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
 par(mfrow=c(3,7))
   plot(seq(1,R,1),M.E[,8], type="l",lwd=1,xlab="Iteration",  ylab=expression(p[12]),main=expression(e[12]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,15], type="l",lwd=1,xlab="Iteration", ylab=expression(p[13]),main=expression(e[13]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,16], type="l",lwd=1,xlab="Iteration", ylab=expression(p[23]),main=expression(e[23]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,22], type="l",lwd=1,xlab="Iteration", ylab=expression(p[14]),main=expression(e[14]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,23], type="l",lwd=1,xlab="Iteration", ylab=expression(p[24]),main=expression(e[24]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,24], type="l",lwd=1,xlab="Iteration", ylab=expression(p[34]),main=expression(e[34]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,28], type="l",lwd=1,xlab="Iteration", ylab=expression(p[15]),main=expression(e[15]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,30], type="l",lwd=1,xlab="Iteration", ylab=expression(p[25]),main=expression(e[25]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,31], type="l",lwd=1,xlab="Iteration", ylab=expression(p[35]),main=expression(e[35]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,32], type="l",lwd=1,xlab="Iteration", ylab=expression(p[45]),main=expression(e[45]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,36], type="l",lwd=1,xlab="Iteration", ylab=expression(p[16]),main=expression(e[16]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,37], type="l",lwd=1,xlab="Iteration", ylab=expression(p[26]),main=expression(e[26]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,38], type="l",lwd=1,xlab="Iteration", ylab=expression(p[36]),main=expression(e[36]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,39], type="l",lwd=1,xlab="Iteration", ylab=expression(p[46]),main=expression(e[46]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,40], type="l",lwd=1,xlab="Iteration", ylab=expression(p[56]),main=expression(e[56]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,43], type="l",lwd=1,xlab="Iteration", ylab=expression(p[17]),main=expression(e[17]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,44], type="l",lwd=1,xlab="Iteration", ylab=expression(p[27]),main=expression(e[27]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,45], type="l",lwd=1,xlab="Iteration", ylab=expression(p[37]),main=expression(e[37]==0),cex.axis=.8)
   plot(seq(1,R,1),M.E[,46], type="l",lwd=1,xlab="Iteration", ylab=expression(p[47]),main=expression(e[47]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,47], type="l",lwd=1,xlab="Iteration", ylab=expression(p[57]),main=expression(e[57]==1),cex.axis=.8)
   plot(seq(1,R,1),M.E[,48], type="l",lwd=1,xlab="Iteration", ylab=expression(p[67]),main=expression(e[67]==1),cex.axis=.8)

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
print("*******************  R E S U L T  *******************")
#****************************************************************

for(cut in c(0.4,0.5,0.6)){
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

#************    Run an example     ********** 

#********************************************* 


DAG(simulate=200,burnin=50,Break=2)






