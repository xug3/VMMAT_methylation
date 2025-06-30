#' Calculate prospective and retrospective P-values for RVMMAT test
#'
#' This function tests a SNPs for a given SNP set for a given estimated null model.
#'
#' @param rvmmat.est The output of function "rvmmat_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants.
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#'
#' @return This function returns a dataframe. The row name is the SNP ID, the first column is the prospective pvalue and the second colum is the restrospective pvalue.
#'
#' @export


vmmat_test <-function(vmmat.est, G)
{
  res<-vmmat.est$Y.res; phi=vmmat.est$phi; V.inv<-vmmat.est$V.inv;X<-vmmat.est$X;N<-nrow(X)
  m<-vmmat.est$m;time<-vmmat.est$time;mu<-vmmat.est$mu;tau<-vmmat.est$tau;cluster.id<-vmmat.est$cluster.id
  snp.names<-colnames(G); family = vmmat.est$family;Rr=vmmat.est$Rr;inciN=vmmat.est$inciN;Td=vmmat.est$Td;Hd=vmmat.est$Hd
  ntime=length(unique(time[,2]))
  
  if(is.vector(G))
  {
   p=1
   }else
  {
    p=dim(G)[2]
  }
  
  P1<-vmmat.est$P1; P2<-vmmat.est$P2;
  
  ResN=c(res)*inciN;  ResNTd=ResN%*%Td
  BTResNTd=matrix(0,nrow=m,ncol=Hd)
  
  IDtable=table(time[,1])
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNTd[ss,]= ResNTd[which(time[,1]==ss),]
    }else{ BTResNTd[ss,]= colSums(as.matrix(ResNTd[which(time[,1]==ss),]))
    }
  }
  
  halfR=Re(expm::sqrtm(Rr)); ResNR=ResN%*%halfR;
  BTResNR=matrix(rep(0,m*ntime),ncol=ntime)
  
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNR[ss,]= ResNR[which(time[,1]==ss),]
    }else{ BTResNR[ss,]= colSums(ResNR[which(time[,1]==ss),])
    }
  }
  BTResNR2colS=colSums(BTResNR^2)
  
  type1result=NULL
  for( j in 1:p){
    if(is.vector(G)) {
      GN=G*inciN
    }else{
      GN=G[,j]*inciN
    }

    GNTd=GN%*%Td; colnames(GNTd)=1:Hd; GNR=GN%*%Rr
    tranG=GN%*%halfR
    ZGNTd=GNTd- P2%*%(P1%*%GNTd); ZGNR=tranG-P2%*%(P1%*%tranG)
    V.pro1<-t(ZGNTd)%*%V.inv%*%ZGNTd;
    V.Rr=t(ZGNR)%*%V.inv%*%ZGNR;
    tran_res<-res
    V.pro1.inv=try(solve(V.pro1))
    if('try-error' %in% class(V.pro1.inv)){
        pSmoothsACAT2=NA
    }else{
        score2=(t(tran_res)%*%GNTd%*%V.pro1.inv%*%t(GNTd)%*%tran_res)[1,1]/(phi^2)
      pSmoothsL2=stats::pchisq(score2,df=Hd,lower.tail=F)
      HalfTSKAT=t(tranG)%*%tran_res
      TSKAT=sum(HalfTSKAT^2)
      lambdaS=eigen(V.Rr,symmetric = TRUE,only.values = TRUE)$values
    
      if(sum(lambdaS)==0) { pSmoothsSKAT2=1
      }else{ pSmoothsSKAT2=generalchisq(lambdaS,TSKAT/(phi^2))}
      Tp=(tan((0.5-pSmoothsL2)*pi)+tan((0.5-pSmoothsSKAT2)*pi))/2
      pSmoothsACAT2=0.5-atan(Tp)/pi
      if(pSmoothsACAT2==0) pSmoothsACAT2=2*pSmoothsL2*pSmoothsSKAT2/(pSmoothsL2+pSmoothsSKAT2)
    }
    
    result=pSmoothsACAT2
    type1result=c(type1result,result)
  }
  names(type1result)=snp.names;
  return(type1result)
  
}



glmm_test <-function(glmm.est, G)
{
    res<-glmm.est$Y.res; V<-glmm.est$V; V.inv<-glmm.est$V.inv;X<-glmm.est$X;N<-nrow(X); 
    m<-glmm.est$m;time<-glmm.est$time;mu<-glmm.est$mu;tau<-glmm.est$tau;cluster.id<-glmm.est$cluster.id;
    snp.names<-colnames(G); family = glmm.est$family;

    P1<-glmm.est$P1; P2<-glmm.est$P2;
    Z<-G-P2%*%(P1%*%G);
    #calculate score;
    n.total<-1;n.rep<-as.numeric(table(time[,1]));tran_res<-rep(0, N);V.pro<-0;
    for(i in 1:m)
    {
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        V.inv.i<-V.inv[[i]];mu.i<-mu[index];G.i<-Z[index,];
        if(family$family == "binomial"){
          if(ni>1){
          delta<-diag(mu.i*(1-mu.i))}else{
          delta<-matrix(mu.i*(1-mu.i))
          }
        }else{
          if(ni>1){
          delta  <- diag(length(mu.i));}else{
          delta <- matrix(1);
          }
        }

        tran_res[index]<-res[index];
        if(ni>1){
          V.pro<-V.pro+t(G.i)%*%V.inv.i%*%G.i;}else{
          V.pro<- V.pro + G.i%*%V.inv.i%*%(G.i)
        }
        phi=glmm.est$phi

    }

    score = c(t(G)%*%tran_res);
    std.pro<-sqrt(diag(V.pro));
    score.pro<-score/(std.pro*phi);
    pval.pro<-pchisq((score.pro)^2, df=1,lower.tail = F)

    result<-cbind(score.pro, pval.pro)
    rownames(result)=snp.names;
    result <- as.data.frame(result)
    return(result)

}


generalchisq=function(lambda,Q){
  muq=sum(lambda)
  sigmaq=sqrt(2*sum(lambda^2))
  s1=sum(lambda^3)/sum(lambda^2)^1.5
  s2=sum(lambda^4)/sum(lambda^2)^2
  if((s1^2)>s2){
    a=1/(s1-sqrt(s1^2-s2))
    delta=s1*a^3-a^2
    l=a^2-2*delta
  }else{
    delta=0
    l=sum(lambda^2)^3/sum(lambda^3)^2
  }
  stats::pchisq((sum(Q)-muq)/sigmaq*sqrt(2*(l+2*delta))+l+delta,df=l, ncp=delta,lower.tail=FALSE)}

