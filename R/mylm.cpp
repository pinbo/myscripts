/*********************
Junli Zhang, July 2015
Modified fastLm for single marker regression
sources for credits:
    * fastLm: https://github.com/RcppCore/RcppArmadillo/blob/master/R/fastLm.R
    * cleanmat: http://stackoverflow.com/questions/13755547/fastest-way-to-drop-rows-with-missing-values
    * fastLm vs lm: http://stackoverflow.com/questions/20034737/rcpparmadillo-fastlm-results-differ-from-rs-lm-what-have-i-done-wrong
    * cbind: http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-February/007139.html
Other useful informations:
    * http://dirk.eddelbuettel.com/code/rcpp.html
    * http://adv-r.had.co.nz/Rcpp.html
    * 
	
To use:
library(Rcpp)
sourceCpp("mylm.cpp")
mylm(single_trait_vector, marker_matrix)
*********************/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// remove rows with missing data from a matrix
// using arma functions
mat cleanmat(NumericMatrix Xr) {
    // get dimensions
    int n = Xr.nrow(),k = Xr.ncol();
    mat X(Xr.begin(), n, k, false );
        // create keep vector
    vec keep = ones<vec>(n);
    for (int j = 0; j < k; j++) 
        for (int i = 0; i < n; i++) 
            if (keep[i] && !is_finite(X(i,j))) keep[i] = 0;
    return X.rows(find(keep==1));
}


// cbind: add a intercept column (all are 1) between y and x
// output for cleanmat to remove rows with missing values
NumericMatrix add1 (NumericVector x, NumericVector y){
	int n = x.size();
	NumericVector z(n, 1.0); // with a default value 1
	NumericMatrix out (x.size(), 3);
	out(_,0) = y; out(_,1)=z; out(_,2)=x;
	return out;
}

// fast linear model using RcppArmadillo functions
// it does not include intercept by default, so I have to cbind(1,x)
List fastLm1(NumericVector yr, NumericMatrix Xr) {
	int n = Xr.nrow(), k = Xr.ncol();
        int df = n - k;
	mat X(Xr.begin(), n, k, false );
	colvec y(yr.begin(), yr.size(), false );
	colvec coef = solve(X, y);
	colvec resid = y - X*coef;
	double sig2 = as_scalar(trans(resid)*resid/(n-k));
	colvec stderrest = sqrt(
		sig2 * diagvec( inv(trans(X)*X)) );
	return List::create(Named("coefficients") = coef,
                Named("stderr") = stderrest,
                Named("df.residual") = df
                );
}


/*
 * work flow:
 * 1) cbind y,1,x using add1 function
 * 2) remove rows with missing values with cleanmat function
 * 3) split y, 1, x again for fastLm1 input
*/




// for multiple markers, but only for one traits input, for multiple traits using version 2 or
// cc <- lapply(all[c(2,3)], function(x) mylm(x,marker3))
// [[Rcpp::export]]
List mylm(NumericVector y, NumericMatrix x){
	//int n = x.size(); //for a list
    int n = x.ncol();
    List out(n);
	for(int i = 0; i < n; ++i) {
        //NumericMatrix m1 = add1(x[i], y);
        NumericMatrix m1 = add1(x(_,i), y);
        NumericMatrix mm = wrap(cleanmat(m1)); //wrap for arma::mat to Rcpp::NumericMatrix
        int m = mm.nrow();
        // seems x.column(i) and x(_,i) can both extract columns
        NumericVector y2 = mm.column(0);
        NumericMatrix x2 (m, 2);
        x2(_,0) = mm.column(1); x2(_,1) = mm.column(2);
        out[i] = fastLm1(y2, x2);
        //out[i] = x2;
	}
	return out;
}

