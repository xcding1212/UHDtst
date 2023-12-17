
TwoWaySampling = function(X, Y, n){
    n1 = nrow(X)
    n2 = nrow(Y)
    if(n>min(max(n1,n2)/2,n1,n2)){
        stop("Invalid sample sizes!")
    }
    if(n1<=n2){
        xs_idx = sample(n1, n)
        ys_idx = sample(n2, 2*n)
        zs_idx = ys_idx[(n+1):(2*n)]
        ys_idx = ys_idx[1:n]
        Xs = X[xs_idx,]
        Ys = Y[ys_idx,]
        Zs = Y[zs_idx,]
    }else{
        xs_idx = sample(n1, 2*n)
        ys_idx = sample(n2, n)
        zs_idx = xs_idx[(n+1):(2*n)]
        xs_idx = xs_idx[1:n]
        Xs = X[xs_idx,]
        Ys = Y[ys_idx,]
        Zs = X[zs_idx,]
    }
    return(list(Xs=Xs, 
                Ys=Ys,
                Zs=Zs))
}


# find eigenvalues of sample covariance matrix under different normalization term
CovEigs = function(X){
    n = nrow(X)
    p = ncol(X)
    X1 = X[]-colMeans(X[])[col(X[])]
    return(eigen(X1%*%t(X1)/sqrt(n*p))$values)
}


K_ = function(x){
    if(abs(x)>=1.05){
        return(0)
    }else if(abs(x)<=1){
        return(1)
    }else{
        return(exp(1/0.05^2-1/(0.05^2-(abs(x)-1)^2)))
    }
}


T_ = function(lambda, gamma, eta0){
    sum = 0
    for(i in 1:length(lambda)){
        sum = sum + (lambda[i]-gamma)/eta0*K_((lambda[i]-gamma)/eta0)
    }
    return(sum) 
}


# check whether our data is of efficient splitting
check_efficient = function(gamma, lambda1, lambda2, epsilon=0.05){
    range1 = abs(lambda1[1]-tail(lambda1,1))
    range2 = abs(lambda2[1]-tail(lambda2,1))
    if(max(abs(gamma-lambda1[1]),abs(gamma-tail(lambda1,1)))>range1-epsilon){
        return(FALSE)
    }else if(max(abs(gamma-lambda2[1]),abs(gamma-tail(lambda2,1)))>range2-epsilon){
        return(FALSE)
    }else{
        return(TRUE)
    }
}


# Put add together, run for once
TwoSampleTest_ = function(X, Y, n, const=0.5, alpha=0.05, epsilon=0.05, thres=NA, mode="test"){
    # split the data
    sample_list = TwoWaySampling(X, Y, n)
    Xs = sample_list$Xs
    Ys = sample_list$Ys
    Zs = sample_list$Zs
    
    # compute eigenvalues and check if splitting is efficient
    eig_xs = CovEigs(Xs)
    eig_ys = CovEigs(Ys)
    eig_zs = CovEigs(Zs)
    gamma = median(eig_zs)
    if(check_efficient(gamma, eig_xs, eig_ys, epsilon)==FALSE){
        return(list(efficient=FALSE,
                    c=1))
    }
    eig_zs = CovEigs(Zs)
    gamma = median(eig_zs)
    eta0 = sd(eig_zs)*const
    
    # compute Tx and Ty
    Tx = T_(eig_xs, gamma, eta0)
    Ty = T_(eig_ys, gamma, eta0)
    if(mode!="test"){
        return(abs(Tx-Ty))   
    }
    
    # compute variance
    if(is.na(thres)){
        thres = 2.6
    }
    
    # rejection rule 1 for rejection, 0 for accept
    if(abs(Tx-Ty)>(thres/qnorm(1-0.05/2)*qnorm(1-alpha/2))){
        return(list(efficient=TRUE,
                    c=1,
                    statistic=(Tx-Ty)/(thres/qnorm(1-0.05/2))))
    }else{
        return(list(efficient=TRUE,
                    c=0,
                    statistic=(Tx-Ty)/(thres/qnorm(1-0.05/2))))
    }
}


#' @export
TwoSampleTest = function(X, Y, n=NA, k=100, const=NA, alpha=0.05, epsilon=0.05, thres=NA, calib=FALSE){
    reject = 0
    df = 0
    statistic = 0
    if(is.na(n)){
        n1 = nrow(X)
        n2 = nrow(Y)
        n = min(max(n1,n2)/2,n1,n2)-5
    }
    if(is.na(thres)){
        if(calib==FALSE){
            thres = 2.6
        }else{
            thres = calibration(n1=nrow(X), n2=nrow(Y), p=ncol(X), n=n, alpha=alpha, const=0.5, iter=100)
        }
    }
    if(is.na(const)){
        const = C_tuning(X, Y, n, alpha=alpha, epsilon=epsilon, thres=thres)$c
    }
    for(i in 1:k){
        result = TwoSampleTest_(X, Y, n, alpha=alpha, const=const, epsilon=epsilon, thres=thres)
        if(result$efficient==TRUE){
            df = df + 1
            statistic = statistic + result$statistic^2
        }
        if(result$c==1){
            reject = reject + 1
        }
    }
    if(df<10){
        pvalue = 0
        statistic = 10000
    }else{
        pvalue = 1 - pchisq(statistic, df)
    }
    # 1 for reject, 0 for accept
    # if(reject>qbinom(1-alpha, k, alpha)){
    #     return(list(decision=1,
    #                 statistic=statistic,
    #                 pvalue=pvalue,
    #                 df=df,
    #                 reject=reject))
    # }else{
    #     return(list(decision=0,
    #                 statistic=statistic,
    #                 pvalue=pvalue,
    #                 df=df,
    #                 reject=reject))
    # }
    if(reject>qbinom(1-alpha, k, alpha)){
        return(list(decision=1,
                    df=df,
                    reject=reject))
    }else{
        return(list(decision=0,
                    df=df,
                    reject=reject))
    }
}


calibration = function(n1, n2, p, n, alpha=0.05, const=0.5, iter=100, K=100){
    values = c()
    for(i in 1:iter){
        X = matrix(rnorm(n1*p, mean=0, sd=1), nrow=n1, ncol=p)
        Y = matrix(rnorm(n2*p, mean=0, sd=1), nrow=n2, ncol=p)
        for(j in 1:K){
            values = c(values, TwoSampleTest_(X, Y, n, alpha=alpha, const=const, mode="calib"))
        }
    }
    return(quantile(values, 1-alpha, names=FALSE))
}



C_tuning = function(X, Y, n, thres=NA, alpha=0.05, epsilon=0.05, K=500){
    Cs = (1:30)/10
    rates = c()
    if(is.na(thres)){
        thres = 2.6
    }
    for(i in 1:length(Cs)){
        all = 0
        rej = 0
        for(j in 1:K){
            result = TwoSampleTest_(X, Y, n, alpha=alpha, const=Cs[i], epsilon=epsilon, thres=thres)
            if(result$efficient){
                all = all + 1
                rej = rej + result$c
            }
        }
        if(all>50){
            rates = c(rates, rej/all)
        }else{
            rates = c(rates, 1)
        }
    }
    stable = find_stable(rates)
    return(list(c=Cs[stable],
                rates=rates))
}

find_stable = function(xs){
    # Here we choose the first drop in variance
    roll_average = zoo::rollmean(xs,3)
    vars = movevar(roll_average)
    for(i in 1:(length(vars)-1)){
        if((vars[i+1]<vars[i])&(roll_average[i])>max(roll_average)/5){
            break
        }
    }
    return(i+2)
}

movevar = function(xs){
    n= length(xs)
    vars = c()
    for(i in 2:n){
        vars = c(vars, var(xs[1:i]))
    }
    return(vars)
}



##################################################
################ Other Methods ###################
##################################################


#' @export
CLX2013 = function(X,Y){
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    p <- dim(X)[2]
    
    W1 <- X[]-colMeans(X[])[col(X[])]
    W2 <- Y[]-colMeans(Y[])[col(Y[])]
    
    S1 <- t(W1) %*% W1 / n1
    S2 <- t(W2) %*% W2 / n2
    
    Theta1 <- Theta2 <- matrix(0,p,p)
    
    for(i in 1:n1){
        Theta1 <- Theta1 + (1/n1) * (  W1[i,] %*% t(W1[i,]) - S1 )^2
    }
    for(i in 1:n2){
        Theta2 <- Theta2 + (1/n2) * (  W2[i,] %*% t(W2[i,]) - S2 )^2
    }
    
    W <- ( S1 - S2 ) / sqrt(Theta1/n1 + Theta2/n2)
    M <- W^2
    M.n <- max(M)
    # qqnorm(W)
    TSvalue <- M.n - 4*log(p) + log(log(p))
    pvalue <- 1 - exp( - 1/sqrt(8*pi) * exp(-TSvalue/2))
    
    output <- list(TSvalue=TSvalue,
                   pvalue=pvalue)
    return(output)
}


#' @export
SY2010 = function(X, Y){
    n1 = nrow(X)-1
    n2 = nrow(Y)-1
    m = ncol(X)
    n = n1 + n2
    
    # try to accelerate the computation
    X = X[]-colMeans(X[])[col(X[])]
    Y = Y[]-colMeans(Y[])[col(Y[])]
    if(m<=n){
        V1 = t(X) %*% X
        V2 = t(Y) %*% Y
    }else{
        V1 = X %*% t(X)
        V2 = Y %*% t(Y)
    }
    V = t(X) %*% X + t(Y) %*% Y
    
    a11 = sum(diag(V1))/(m*n1)
    a12 = sum(diag(V2))/(m*n2)
    # since tr(V1V1^T)=||V1||_F^2
    a21 = (norm(V1, type="F")^2-sum(diag(V1))^2/n1)/(m*(n1-1)*(n1+2))
    a22 = (norm(V2, type="F")^2-sum(diag(V2))^2/n2)/(m*(n2-1)*(n2+2))
    
    # rank(V)<= min(2*m,n1+n2+2)
    r = min(m, n1+n2+2) # largest rank of V
    eigV = suppressWarnings(rARPACK::eigs(V, k=r)$values) 
    a1 = sum(eigV)/(n*m)
    a2 = (sum(eigV^2)-sum(diag(V))^2/n)/(m*(n-1)*(n+2))
    a3 = (sum(eigV^3)/m-3*n*(n+1)*m*a2*a1-n*m^2*a1^3)/(n*(n^2+3*n+4))
    c0 = n*(n^3+6*n^2+21*n+18)
    c1 = 2*n*(2*n^2+6*n+9)
    c2 = 2*n*(3*n+2)
    c3 = n*(2*n^2+5*n+7)
    a4 = (sum(eigV^4)/m-m*c1*a1-m^2*c2*a1^2*a2-m*c3*a2^2-n*m^3*a1^4)/c0
    
    gamma1 = a21/a11^2
    gamma2 = a22/a12^2
    xi1_2 = 4/n1^2*(a2^2/a1^4+2*n1/m*(a2^3/a1^6-2*a2*a3/a1^5+a4/a1^4))
    xi2_2 = 4/n2^2*(a2^2/a1^4+2*n2/m*(a2^3/a1^6-2*a2*a3/a1^5+a4/a1^4))
    
    Q2 = (gamma1-gamma2)^2/(xi1_2+xi2_2)
    pvalue = pchisq(Q2, df=1, lower.tail=FALSE)
    return(list(Q2=Q2,
                pvalue=pvalue))
}


#' @export
LC2012 = function(X, Y){
    n1 = nrow(X)
    n2 = nrow(Y)
    result = equalCovs::equalCovs(X,Y,n1,n2)
    return(list(statistic=result[1],
                pvalue=result[2]))
}



#' @export
HC2018 = function(X, Y, N=floor(ncol(X)^(0.7)), alpha=0.05){
    double_sum = function(X1,X2){
        colSums(X1)*colSums(X2)-colSums(X1*X2)
    }
    triple_sum = function(X1,X2,X3){
        double_sum(X1,X2)*colSums(X3)-double_sum(X1*X3,X2)-double_sum(X1,X2*X3)
    }
    quad_sum = function(X1,X2,X3,X4){
        triple_sum(X1,X2,X3)*colSums(X4)-triple_sum(X1*X4,X2,X3)-triple_sum(X1,X2*X4,X3)-triple_sum(X1,X2,X3*X4)
    }
    Di = function(X, q){
        # We here optimize the computation of U-statistics, we try to avoid for-loops
        n = nrow(X)
        p = ncol(X)
        X1 = X[,1:(p-q)]
        X2 = X[,(q+1):p]
        # use some product to replace multiple summation
        # we use inclusion-exclusion principle to deal with this
        D_1 = sum(double_sum(X1*X2,X1*X2))
        D_2 = sum(triple_sum(X1,X2,X1*X2))
        D_3 = sum(quad_sum(X1,X2,X1,X2))
        1/(n*(n-1))*D_1-2/(n*(n-1)*(n-1))*D_2+1/(n*(n-1)*(n-2)*(n-3))*D_3
    }
    Dc = function(X1, X2, q){
        n1 = nrow(X1)
        n2 = nrow(X2)
        p = ncol(X1)
        X11 = X1[,1:(p-q)]
        X12 = X1[,(q+1):p]
        X21 = X2[,1:(p-q)]
        X22 = X2[,(q+1):p]
        Dc_1 = sum(colSums(X11*X12)*colSums(X21*X22))
        Dc_2 = sum(double_sum(X11,X12)*colSums(X21*X22))
        Dc_3 = sum(colSums(X11*X12)*double_sum(X21,X22))
        Dc_4 = sum(double_sum(X11,X12)*double_sum(X21,X22))
        Dc_1/(n1*n2)-Dc_2/(n1*(n1-1)*n2)-Dc_3/(n1*n2*(n2-1))+Dc_4/(n1*(n1-1)*n2*(n2-1))
    }
    Sq = function(X1,X2,q){
        Di(X1,q)+Di(X2,q)-2*Dc(X1,X2,q)
    }
    Ri = function(X, q){
        n = nrow(X)
        p = ncol(X)
        X = X[]-colMeans(X[])[col(X[])]
        X1 = X[,1:(p-q)]
        X2 = X[,(q+1):p]
        Y = X1*X2
        Y = Y-colSums(Y[])[col(Y[])]/(n-1)
        YYt2 = (Y%*%t(Y))^2
        (sum(YYt2)-sum(diag(YYt2)))/(n*(n-1))
    }
    Rc = function(X1,X2,q){
        n1 = nrow(X1)
        p = ncol(X1)
        X1 = X1[]-colMeans(X1[])[col(X1[])]
        X11 = X1[,1:(p-q)]
        X12 = X1[,(q+1):p]
        Y1 = X11*X12
        Y1 = Y1-colSums(Y1[])[col(Y1[])]/(n1-1)
        
        n2 = nrow(X2)
        X2 = X2[]-colMeans(X2[])[col(X2[])]
        X21 = X2[,1:(p-q)]
        X22 = X2[,(q+1):p]
        Y2 = X21*X22
        Y2 = Y2-colSums(Y2[])[col(Y2[])]/(n2-1)
        sum((Y1%*%t(Y2))^2)/(n1*n2)
    }
    V2 = function(X1,X2,q){
        n1 = nrow(X1)
        n2 = nrow(X2)
        Ri(X1,q)*2/(n1*(n1-1))+Ri(X2,q)*2/(n2*(n2-1))+Rc(X1,X2,q)*4/(n1*n2)
    }
    one_super = function(X1,X2,q){
        chi = Sq(X1,X2,q)^2/V2(X1,X2,q)
        pchisq(chi,1,lower.tail=FALSE)
    }
    pvalues = c()
    for(i in 0:N){
        pvalues = c(pvalues, one_super(X,Y,i))
    }
    test = mutoss::adaptiveSTS(pvalues, alpha=alpha, lambda=0.5, silent=TRUE)
    return(list(reject=sum(test$rejected),
                pvalues=pvalues,
                N=N))
}
