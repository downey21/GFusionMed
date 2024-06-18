#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericMatrix t_cpp(NumericMatrix mat) {
    int rows = mat.nrow();
    int cols = mat.ncol();

    NumericMatrix result(cols, rows);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(j, i) = mat(i, j);
        }
    }

    return result;
}

// [[Rcpp::export]]
NumericMatrix prod_cpp(NumericMatrix mat1, NumericMatrix mat2) {
    int nrow1 = mat1.nrow();
    int ncol1 = mat1.ncol();
    int ncol2 = mat2.ncol();

    NumericMatrix result(nrow1, ncol2);

    for (int i = 0; i < nrow1; ++i) {
        for (int j = 0; j < ncol2; ++j) {
            double sum = 0.0;
            for (int k = 0; k < ncol1; ++k) {
                sum += mat1(i, k) * mat2(k, j);
            }
            result(i, j) = sum;
        }
    }

    return result;
}

// [[Rcpp::export]]
NumericMatrix t_prod_cpp(NumericMatrix mat1, NumericMatrix mat2) {
    NumericMatrix transposedMat1 = t_cpp(mat1);
    NumericMatrix result = prod_cpp(transposedMat1, mat2);

    return result;
}

// [[Rcpp::export]]
NumericMatrix prod_t_cpp(NumericMatrix mat1, NumericMatrix mat2) {
    NumericMatrix transposedMat2 = t_cpp(mat2);
    NumericMatrix result = prod_cpp(mat1, transposedMat2);

    return result;
}

// [[Rcpp::export]]
NumericMatrix inverse_cpp(NumericMatrix mat) {
    const Map<MatrixXd> matEigen(as<Map<MatrixXd> >(mat));
    MatrixXd invMatEigen = matEigen.inverse();
    return wrap(invMatEigen);
}