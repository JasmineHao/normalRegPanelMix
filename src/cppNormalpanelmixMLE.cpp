#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::export]]

List cppregpanelmixMLE(NumericVector ys,
                            NumericMatrix xs,
                            NumericMatrix zs,
                           NumericVector alpha0s,
                           NumericVector mu0s,
                           NumericVector sigma0s,
                           NumericMatrix beta0s,
                           NumericVector gamma0s,
                           int m,
                           int q,
                           int p,
                           int t){
  int nt = ys.size();
  int n  = nt / t;
  vec alpha(alpha0s.begin(), alpha0s.size(), false);
  vec y(ys.begin(), ys.size(), false);
  mat x(xs.begin(), xs.nrow(),xs.ncol(),false);
  mat z(zs.begin(), zs.nrow(), zs.ncol(), false);
  vec mu(mu0s.begin(), mu0s.size(), false);
  vec sigma(sigma0s.begin(), sigma0s.size(), false);
  vec gamma(gamma0s.begin(), gamma0s.size(),false);
  mat beta(beta0s.begin(),beta0s.nrow(),beta0s.ncol(),false);
  int q1 = q+1;
  mat mubeta(q1,m);
  // mat mubeta(q1,m);
  vec r(m), r_t(m),l_j(m); //NEED this r_i since need to sum the squared residuals
  mat w(m,nt); //Weight matrix, note weight will be similar for same i
  vec ytilde(nt);
  vec wtilde(nt);
  mat ztilde(nt,p);
  mat x1(nt,q1);

  double minr, sum_l_j;
  double ll = 0; // force initilization
  ytilde = y;
  ll = - (double)nt * M_LN_SQRT_2PI; /* n/2 times log(2pi) */
  for (int j=0; j < m; ++j){
    mubeta(0,j) = mu(j);
    for (int i=0; i<q; ++i){
      mubeta(i+1,j) = beta(j,i);
    }
  }
  x1.col(0) = ones(nt);
  for (int i = 1; i < q1; ++i){
    x1.col(i) = x.col(i-1);
  }
  for (int nn = 0; nn < n; nn++) {
    r.fill(0.0);
    for (int tt = 0; tt < t; tt++){
      /* standardized squared residual */
      /* posterior for i */
      // Rcout << trans(x1.row(nn*t + tt)*mubeta) << endl;
      r_t = (1.0/sigma) % (ytilde(nn*t+tt) - trans(x1.row(nn*t + tt)*mubeta));
      r_t = 0.5 * (r_t % r_t); /* This is faster than r = pow( r, 2.0 ) */
      r += r_t; //sum the residual for each time period

    }
    r = r + t*log(sigma);



    /* normalizing with minr avoids the problem of dividing by zero */
    minr = min(r);
    // Rcout << exp(r - minr) << endl;
    l_j =  alpha % exp( minr-r );
    sum_l_j = sum( l_j );
    for (int tt= 0; tt<t;tt++){
      w.col(nn*t + tt) = l_j/sum_l_j; /* w(j,i) = alp_j*l_j / sum_j (alp_j*l_j) */
    }
    /* loglikelihood*/
    ll +=  log(sum_l_j) - minr; /* subtract back minr */
  } /* end for (i=0; i<n; i++) loop */

  return List::create(Named("likelihood") = ll,
                            Named("Dummy") = 1.0);

}
