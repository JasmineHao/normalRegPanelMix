#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double SINGULAR_EPS = 10e-10; // criteria for matrix singularity

// [[Rcpp::export]]
List cppRegPanelmixPMLE(NumericMatrix bs,
                           NumericVector ys,
                           NumericMatrix xs,
                           NumericMatrix zs,
                           NumericVector mu0s,
                           NumericVector sigma0s,
                           int m,
                           int p,
                           int t, //time periods,crucial for mix
                           double an,
                           int maxit = 2000,
                           int ninits = 10,
                           double tol = 1e-8,
                           double tau = 0.5,
                           int h = 0,
                           int k = 0,
                           int update_alpha=1) {
  int nt = ys.size();
  int n  = nt / t;
  int q = xs.ncol();
  int q1 = q + 1;
  arma::mat b(bs.begin(), bs.nrow(), bs.ncol(), false);
  arma::vec y(ys.begin(), ys.size(), false);
  arma::mat x(xs.begin(), xs.nrow(), xs.ncol(), false);
  arma::mat z(zs.begin(), zs.nrow(), zs.ncol(), false);
  arma::vec mu0(mu0s.begin(), mu0s.size(), false);
  arma::vec sigma0(sigma0s.begin(), sigma0s.size(), false);
  arma::vec b_jn(bs.nrow());
  arma::vec lb(m),ub(m);
  arma::vec alpha(m), mu(m), sigma(m), alp_sig(m);
  arma::mat mubeta(q1,m);
  arma::vec r(m), r_t(m),l_j(m); //NEED this r_i since need to sum the squared residuals
  arma::mat w(m,nt); //Weight matrix, note weight will be similar for same i
  arma::mat post(m*nt,ninits);
  arma::vec notcg(ninits), penloglikset(ninits), loglikset(ninits);
  arma::vec gamma(p);
  arma::vec ytilde(nt);
  arma::mat xtilde(nt,q1);
  arma::vec wtilde(nt);
  arma::mat ztilde(nt,p);
  arma::mat zz(p,p);
  arma::mat ze(p,1);
  int sing, emit;
  double oldpenloglik, s0j, diff, minr, w_j, sum_l_j, ssr_j, alphah, tauhat;
  arma::mat x1(nt,q1);
  double ll = 0; // force initilization
  double penloglik = 0; // force initialization

  x1.col(0) = arma::ones(nt);
  for (int i=1; i < q1; ++i) {
    x1.col(i) = x.col(i-1);
  }
  /* Lower and upper bound for mu */
  if (k==1) {  // If k==1, compute upper and lower bounds
    mu0(0) = R_NegInf;
    mu0(m) = R_PosInf;
    for (int j=0; j<h; j++) {
      lb(j) = (mu0(j)+mu0(j+1))/2.0;
      ub(j) = (mu0(j+1)+mu0(j+2))/2.0;
    }
    for (int j=h; j<m; j++) {
      lb(j) = (mu0(j-1)+mu0(j))/2.0;
      ub(j) = (mu0(j)+mu0(j+1))/2.0;
    }
  }

  /* iteration over ninits initial values of b */
  for (int jn=0; jn<ninits; jn++) {

    /* initialize EM iteration */
    b_jn = b.col(jn);
    for (int j=0; j < m; ++j){
      alpha(j) = b_jn(j);
      for (int i=0; i<q1; ++i){
        mubeta(i,j) = b_jn(m+q1*j+i);
      }
      sigma(j) = b_jn((q1+1)*m+j);
    }
    if (p>0) {
      for (int j=0; j < p; j++){
        gamma(j) = b_jn((q1+2)*m+j);
      }
    }
    oldpenloglik = R_NegInf;
    emit = 0;
    diff = 1.0;
    sing = 0;

    /* EM loop begins */

    for (int iter = 0; iter < maxit; iter++) {
      ll = - (double)nt * M_LN_SQRT_2PI; /* nt/2 times log(2pi) */
      // alp_sig = alpha/pow(sigma,t);

    if (p==0) {
      ytilde = y;
    } else {
      ytilde = y - z*gamma;
    }

    for (int nn = 0; nn < n; nn++) {
      r.fill(0);
      for (int tt = 0; tt < t; tt++){
        /* standardized squared residual */
        r_t = (1.0/sigma) % (ytilde(nn*t + tt) - trans(x1.row(nn*t + tt)*mubeta));
        r_t = 0.5 * (r_t % r_t); /* This is faster than r = pow( r, 2.0 ) */
        r = r + r_t; //sum the residual for each time period
      }
      r = r + t*log(sigma);
      /* normalizing with minr avoids the problem of dividing by zero */
      minr = min(r);
      l_j =  alpha % exp( minr-r );
      sum_l_j = sum( l_j );
      for (int tt= 0; tt<t;tt++){
        w.col(nn*t + tt) = l_j/sum_l_j; /* w(j,i) = alp_j*l_j / sum_j (alp_j*l_j) */
      }
      /* loglikelihood*/
      ll +=  log(sum_l_j) - minr; /* subtract back minr */
    } /* end for (i=0; i<n; i++) loop */

      /* Compute the penalized loglik. Note that penalized loglik uses old (not updated) sigma */
      penloglik = ll + log(2.0) + fmin(log(tau),log(1-tau));
    for (int j=0; j<m; j++) {
      s0j = sigma0(j)/sigma(j);
      penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
    }
    diff = penloglik - oldpenloglik;
    oldpenloglik = penloglik;

    /* Normal exit */
    if (diff < tol ){
      break;
    }
    emit++;

    /* update alpha, mu, and sigma */
    for (int j = 0; j < m; j++) {
      w_j = sum( w.row(j) ); /* w_j(j) = sum_i w(i,j) */
      if (update_alpha == 1){
        alpha(j) = w_j / nt;
      }      
      wtilde = trans(w.row(j));
      for (int ii = 0; ii < q1; ii++) {
        xtilde.col(ii) = wtilde % x1.col(ii);
      }
      /* Check for singularity*/
      arma::mat design_matrix = trans(xtilde) * x1;
      if (rcond(design_matrix) < SINGULAR_EPS)
      {
        sing = 1;
        break;
      }
      mubeta.col(j) = solve( design_matrix , trans(xtilde) * ytilde );

      ssr_j = sum( trans(w.row(j)) % pow(  ytilde - x1*mubeta.col(j), 2 ) );
      sigma(j) = sqrt( (ssr_j + 2.0*an*sigma0(j)*sigma0(j))  / (w_j + 2.0*an) );
      sigma(j) = fmax(sigma(j),0.05*sigma0(j)); //Why changed from 0.01 to 0.05
      /* If k ==1, impose lower and upper bound */
      if (k==1) {
        mubeta(0,j) = fmin( fmax(mubeta(0,j),lb(j)), ub(j));
      }
    }

    /* for PMLE, we set k=0 (default value) */
    /* for EM test, we start from k=1       */
    /*   if k==1, we don't update tau       */
    /*   if k>1, we update tau              */
    if (k==1){
      alphah = (alpha(h-1)+alpha(h));
      alpha(h-1) = alphah*tau;
      alpha(h) = alphah*(1-tau);
    } else if (k>1 && update_alpha == 1) {
      alphah = (alpha(h-1)+alpha(h));
      tauhat = alpha(h-1)/(alpha(h-1)+alpha(h));
      if(tauhat <= 0.5) {
        tau = fmin((alpha(h-1)*n + 1.0)/(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      } else {
        tau = fmax(alpha(h-1)*n /(alpha(h-1)*n + alpha(h)*n + 1.0), 0.5);
      }
      alpha(h-1) = alphah*tau;
      alpha(h) = alphah*(1-tau);
    }

    if (p>0) { /* update gamma */
      zz.zeros();
      ze.zeros();
      for (int j = 0; j < m; j++) {
        wtilde = trans(w.row(j)) ;
        for (int ii = 0; ii < p; ii++) {
          ztilde.col(ii) = wtilde % z.col(ii);
        }
        zz = zz + ( trans(ztilde) * z ) / (sigma(j)*sigma(j));
        ze = ze + ( trans(ztilde) * (y - x1*mubeta.col(j)) ) / (sigma(j)*sigma(j));
      }
      // sanity check before solving an inverse matrix;
      // if it is likely singular, leave it as is.
      if (cond(zz) < SINGULAR_EPS)
      {
        sing = 1;
        break;
      }
      else
        gamma = solve(zz,ze);
    }

    /* Check singularity */
    for (int j=0; j<m; j++) {
      if (alpha(j) < 1e-8 || std::isnan(alpha(j)) || sigma(j) < 1e-8){
        sing = 1;
      }
    }

    /* Exit from the loop if singular */
    if (sing) {
      notcg(jn) = 1;
      break;
    }

    }/* EM loop ends */
    /*print to see if the weight is correct */

    penloglikset(jn) = penloglik;
    loglikset(jn) = ll;
    for (int j=0; j < m; j++){
      b_jn(j) = alpha(j);
      for (int i=0; i<q1; ++i){
        b_jn(m+q1*j+i) = mubeta(i,j);
      }
      b_jn((q1+1)*m+j) = sigma(j);
    }
    if (p>0) {
      for (int j=0; j < p; j++){
        b_jn(3*m+j) = gamma(j);
      }
    }
    b.col(jn) = b_jn; /* b is updated */
    post.col(jn) = vectorise(trans(w));

  } /* end for (jn=0; jn<ninits; jn++) loop */
  return Rcpp::List::create(Named("penloglikset") = wrap(penloglikset),
                            Named("loglikset") = wrap(loglikset),
                            Named("notcg") = wrap(notcg),
                            Named("post") = wrap(post)
  );
}
