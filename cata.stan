/*
  *Simple normal regression example
*/

functions {
  real catapost_lpdf(vector bb, matrix X, matrix Xstar,vector Y, vector Ystar, real tau, int n, int M, int p) {
    vector[n] linpred;
    vector[M] linpredstar;
    real lpf;
    
    linpred=X*bb;
    linpredstar=Xstar*bb;
    lpf=sum(Y .* linpred - log(1+exp(linpred)) ) + tau/M *sum(Ystar .* linpredstar - log(1+exp(linpredstar)) );
    return lpf;
  }
}

 data {
    int n; //the number of observations
    int M; //the size of the pseudo sample
    int p; //the number of columns in the model matrix
    vector[n] Y; //the response
    vector[M] Ystar; //the pseudoresponse
    matrix[n,p] X; //the model matrix
    matrix[M,p] Xstar; //the pseudo matrix
    real tau;  // prior weight
}


parameters {
  vector[p] beta; //the regression parameters
}



model {  
 beta ~ catapost(X, Xstar, Y, Ystar, tau, n, M, p);//binomial(linpred,1);
}
