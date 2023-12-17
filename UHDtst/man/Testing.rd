\name{Testing}
\alias{TwoSampleTest}
\alias{SY2010}
\alias{LC2012}
\alias{CLX2013}
\alias{HC2018}
\title{Ultra-high Dimensional Two Sample Covariance Test}
\description{
    Test of equal covariance matrices of two given samples X and Y
}
\usage{
TwoSampleTest(X, Y, n=NA, k=100, const=NA, alpha=0.05, epsilon=0.05, thres=NA, 
              calib=FALSE)
SY2010(X, Y)
LC2012(X, Y)
CLX2013(X, Y)
HC2018(X, Y, N=floor(ncol(X)^(0.7)), alpha=0.05)
}
\arguments{
    \item{X, Y}{input samples in the form n by p where p is the dimension.}
    \item{n}{numbers of samples in bootstrapping. If it is not given, it will be computed automatically.}
    \item{k}{numbers of bootstrappings.}
    \item{const}{constant for bandwidth. If it is not given, it will be tuned.}
    \item{alpha}{size of our test.}
    \item{epsilon}{constant used to detect efficient splitting.}
    \item{thres}{some tuned constant.}
    \item{calib}{whether to use calibration to tune thres.}
    \item{N}{numbers of super-diagonals to use.}
}
\value{
    \code{TwoSampleTest} tests equal covariance matrices of two samples using DHW2023,

    \code{SY2010} tests equal covariance matrices of two samples using SY2010,

    \code{LC2012} tests equal covariance matrices of two samples using LC2012,
    
    \code{CLX2013} tests equal covariance matrices of two samples using CLX2013,
    
    \code{HC2018} tests equal covariance matrices of two samples using HC2018.
}
\examples{
n1 = 100
n2 = 100
p = 200
X = matrix(rnorm(n1*p), ncol=p)
Y = matrix(rnorm(n2*p), ncol=p)
TwoSampleTest(X, Y, const=0.5)$pvalue
SY2010(X, Y)$pvalue
LC2012(X, Y)$pvalue
CLX2013(X, Y)$pvalue
HC2018(X, Y, alpha=0.05)$reject
}

\references{
    [1] Srivastava, M. S., & Yanagihara, H. (2010). Testing the equality of several covariance matrices with fewer observations than the dimension. Journal of Multivariate Analysis, 101(6), 1319-1329.
    
    [2] Li, J., & Chen, S. X. (2012). TWO SAMPLE TESTS FOR HIGH-DIMENSIONAL COVARIANCE MATRICES. The Annals of Statistics, 40(2), 908–940.
    
    [3] Cai, T., Liu, W., & Xia, Y. (2013). Two-sample covariance matrix testing and support recovery in high-dimensional and sparse settings. Journal of the American Statistical Association, 108(501), 265-277.
    
    [4] He, J., & Chen, S. X. (2018). High-dimensional two-sample covariance matrix testing via super-diagonals. Statistica Sinica, 28(4), 2671-2696.
    
    [5] Ding, X. C. & Hu, Y. C., & Wang, Z. G. (2023). Two sample test for covariance matrices in ultra-high dimension. In progress.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{test}
