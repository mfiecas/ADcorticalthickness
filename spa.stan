functions{
  real flogitg(real alp, real z)
  {
    vector[4] out;
    if (alp==1) return log(exp(z)-1);
    for(i in 1:4) out[i]=(alp-1)*log(z/4*(i-0.5))+z/4*(i-0.5);
    return log_sum_exp(out)-log(4)+log(z);
  }
  
  real nngp_lpdf(vector bet, real phi, real sigmasq, real tausq, int[,] NNind, matrix NNdist, matrix NNdistM, int V, int M)
  {
    vector[V] D; //scaled by sigmasq
    vector[V] U;
    real kappa_p_1=tausq/sigmasq+1;
    int l;
    int h;
    
    D[1]=kappa_p_1;
    U[1]=bet[1];
    for (v in 2:V)
    {
      //l=(v<(M+1)) ? (v-1):M;
      matrix[(v<(M+1)) ? (v-1):M,(v<(M+1)) ? (v-1):M] vNNdistM; //C(N,N)+tausq*I scaled by sigmasq
      matrix[(v<(M+1)) ? (v-1):M,(v<(M+1)) ? (v-1):M] vNNCholL; //cholL of vNNdistM
      vector[(v<(M+1)) ? (v-1):M] vNNcorr; //C(v,N) scaled by sigmasq
      vector[(v<(M+1)) ? (v-1):M] v1;
      row_vector[(v<(M+1)) ? (v-1):M] v2;
      l=(v<(M+1)) ? (v-1):M; //dim/number of nearest neighbors
      
      for(j in 1:l) vNNdistM[j,j]=kappa_p_1;
      if(l>1)
      {
        h=1;
        for(j in 2:l)
        {
          for(k in 1:(j-1)){
            vNNdistM[j,k]=exp(-phi*NNdistM[v-1,h]);
            vNNdistM[k,j]=vNNdistM[j,k];
            h=h+1;
          }
        }
      }
      
      vNNCholL=cholesky_decompose(vNNdistM);
      vNNcorr=to_vector(exp(-phi*NNdist[v-1,1:l]));
      v1=mdivide_left_tri_low(vNNCholL,vNNcorr);
      D[v]=kappa_p_1-dot_self(v1);
      v2 = mdivide_right_tri_low(v1',vNNCholL);
      U[v]=bet[v]-v2*bet[NNind[v-1,1:l]];
    }
    return -0.5*(1/sigmasq*dot_product(U,(U ./ D))+sum(log(D))+V*log(sigmasq));
  }
  
  real surv(vector y, int[] del, vector bet, vector gam, vector alp, matrix Z, matrix x0, matrix k, vector zdel, real logydel, vector x0kydel, int n, int P, int V)
  {
    vector[V] Q0;
    vector[V] Qdel;
    vector[n] Zgam;
    vector[n] x0b;
    vector[n] kb;
    vector[n] logitg;
    
    for(v in 1:V)
    {
      Qdel[v]=(gam[(v-1)*(P+1)+1]+log(alp[v]))*sum(del)+alp[v]*logydel+dot_product(zdel,gam[((v-1)*(P+1)+2):((P+1)*v)])+x0kydel[v]*bet[v];
      
      Zgam=Z*gam[((v-1)*(P+1)+2):((P+1)*v)]+gam[(v-1)*(P+1)+1];
      x0b=x0[,v]*bet[v];
      kb=k[,v]*bet[v];
      for(i in 1:n)
      {
        if (kb[i]==0)
          logitg[i]=alp[v]*log(y[i])-log(alp[v]);
        else if (kb[i]<0)
          logitg[i]=-alp[v]*log(-kb[i])+lgamma(alp[v])+log(gamma_p(alp[v],-kb[i]*y[i]));
        else
          logitg[i]=-alp[v]*log(kb[i])+flogitg(alp[v],kb[i]*y[i]);
      }
      Q0[v]=exp(log_sum_exp(Zgam+x0b+logitg));
    }
    return sum(Qdel)-dot_product(Q0,alp);
  }
}

data
{
  int<lower=1> n; //subjects
  int<lower=0> P; //covariates
  int<lower=1> V; //vertices
  int<lower=1> M; //nearest neighbors
  
  matrix[n,P] Z;
  matrix[n,V] x0;
  matrix[n,V] k;
  vector[n] y;
  int del[n]; //conversion indicator
  
  int NNind[V-1,M];
  matrix[V-1,M] NNdist;
  matrix[V-1,(M*(M-1)/2)] NNdistM;
  
  real<lower=0> sgam;
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> theta1;
  real<lower=0> theta2;
  real<lower=0> ls;
  real<lower=0> lt;
  real<lower=0> ks;
  real<lower=0> kt;
}

transformed data
{
  vector[n] vecdel;
  vector[P] zdel;
  real logydel;
  vector[V] x0kydel;
  vecdel=to_vector(del);
  zdel=Z'*vecdel;
  logydel=dot_product(log(y),vecdel);
  x0kydel=x0'*vecdel+k'*(y .* vecdel);
}

parameters
{
  vector[V] bet;
  vector[(P+1)*V] gam;
  vector<lower=0>[V] alp;
  
  real<lower=theta1,upper=theta2> phi;
  real<lower=0> sigmasq;
  real<lower=0> tausq;
}

model
{
  phi ~ uniform(theta1,theta2);
  sigmasq ~ inv_gamma(ls,ks);
  tausq ~ inv_gamma(lt,kt);
  alp ~ gamma(a,b);
  gam ~ normal(0,sgam);
  bet ~ nngp(phi,sigmasq,tausq,NNind,NNdist,NNdistM,V,M);
  target += surv(y,del,bet,gam,alp,Z,x0,k,zdel,logydel,x0kydel,n,P,V);
}
