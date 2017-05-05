
#' Fits a quantile regression SVM based on the Pinball Loss
#'
#' @param x An n X m matrix containing the predictors (n= number of observatiosn, m = number of predictors) .
#' @param y The Response onto which the qrsvm shall be fitted
#' @param kernel a string giving the type of kernels from package kernlab to use f.e. "rbfdot" for Radial Basis Function Kernel. All Kernels except "stringdot" supported.
#' @param tau The Quantile that shall be estimated. 0<=tau<=1
#' @param sigma A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param degree A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param scale A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param offset A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param order A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand. Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to other SVM implementations.
#' @return An object of class "qrsvm"
#' @references "Nonparametric Quantile Regression" by I.Takeuchi, Q.V. Le, T. Sears, A.J. Smola (2004)
#' @import kernlab
#' @import Matrix
#' @import quadprog
#' @examples
#' n<-300
#' x<-seq(-2,2,length.out=n)
#' x<-as.matrix(x)
#' y<-rnorm(n)*abs(x^2+1)
#'
#' mod1<-qrsvm(x,y)
#' fit1<-mod1$fitted
#'
#' mod2<-qrsvm(x,y, tau=0.05)
#' fit2<-mod2$fitted
#'
#' plot(x,y)
#' lines(x=x, y=fit1, col="red")
#' lines(x=x, y=fit2, col="red")
#' @export
qrsvm <- function(x, y, kernel = "rbfdot", cost = 1,    tau = 0.95,
                  sigma = 5, degree = 2, scale = 1, offset = 1, order = 1) {




    library(kernlab)
    library(quadprog)
    library(Matrix)

    if (kernel == "rbfdot") {
        kern <- rbfdot(sigma = sigma)
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "polydot") {
        kern <- polydot(degree = degree, scale = scale, offset = offset)
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "vanilladot") {
        kern <- vanilladot()
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "tanhdot") {
        kern <- tanhdot(scale = scale, offset = offset)
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "laplacedot") {
        kern <- laplacedot(sigma = sigma)
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "besseldot") {
        kern <- besseldot(sigma = sigma, order = order, degree = degree)
        kernmat <- kernelMatrix(kern, x)
    }
    if (kernel == "anovadot") {
        kern <- anovadot(sigma = sigma, degree = degree)
        kernmat <- kernelMatrix(kern, x)
    }
    if (nrow(kernmat) < nrow(x)) {
        print("Kernelmatrix not valid! Check if valid kernel type stated!")
    }

    Amat <- rbind(rep(1, nrow(x)), diag(x = 1, nrow = nrow(x)),
        diag(x = -1, nrow = nrow(x)))
    pdcalc <- nearPD(kernmat)
    pdmat <- pdcalc$mat
    Dmat <- pdmat
    dvec <- y
    b0 <- c(0, rep((cost * (tau - 1)), nrow(x)), rep(-cost *
        tau, nrow(x)))
    erg <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
        bvec = b0, meq = 1, factorized = FALSE)
    alpha <- erg$solution
      f <- alpha %*% kernmat
    offshift <- which.min((round(alpha, 3) -(cost * tau))^2 + (round(alpha, 3) - (cost * (tau - 1)))^2)
    b <- y[offshift] - f[offshift]
    fnew <- alpha %*% kernmat + b

    model <- list(alpha = alpha, xtrain = x, kernel = kern,
        sigma = sigma, cost = cost, b0 = b, fitted = as.numeric(fnew),
        tau = tau, scale = scale, offset = offset, order = order,
        kernstring = kernel, y = y)
    class(model) <- "qrsvm"
    return(model)
}







