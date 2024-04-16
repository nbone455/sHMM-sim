#  When the true model is a threshold model, the corHMM models are less powerful as they add too many parameters. 
# Heatmap 3: over corHMM model with hidden states (2,4,6,8) when true model is a threshold model. 
# Other axis: Sample size or Ancestral state? Allows direct comparison to heatmap 1. 

## This is why you don't make copy and paste scripts 
## probably better to have all outputs saved in a list of lists: 
## runs 
## (corhmm runs(2,4,6,8,10,16)) (thresh runs(2,4,6,8,10,16))

# load packages and functions
library(ape)
library(treeplyr)
library(textshape)
library(vctrs)
library(stringr)
library(phytools)
require(dplyr)
require(corHMM)
library(geiger)
library(ggplot2)
library(gplots)

# how many simulations
nsim <- 50

#old matrix maker functions from corhmm
rate.mat.maker <- function (rate.cat, hrm = TRUE, ntraits = NULL, nstates = NULL, model = c("ER", "SYM", "ARD")){
  if (hrm == TRUE) {
    k = 2
    mat1 <- matrix(NA, k * rate.cat, k * rate.cat)
    mat2 <- matrix(NA, k * rate.cat, k * rate.cat)
    vec.tmp1 <- rep(c(0, 1), rate.cat)
    vec.tmp2 <- rep(1:rate.cat, rep(2, rate.cat)) - 1
    for (i in 1:(k * rate.cat)) {
      mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
      mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
    }
    matFINAL <- mat1 + mat2
    rate.mat.index <- matrix(NA, k * rate.cat, k * rate.cat)
    np <- k + (rate.cat - 1) * 6
    index <- matFINAL == 1
    rate.mat.index[index] <- 1:np
    if (rate.cat == 1) {
      rownames(rate.mat.index) <- c("(0)", "(1)")
      colnames(rate.mat.index) <- c("(0)", "(1)")
    }
    if (rate.cat == 2) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)")
    }
    if (rate.cat == 3) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)")
    }
    if (rate.cat == 4) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)")
    }
    if (rate.cat == 5) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)")
    }
    if (rate.cat == 6) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)")
    }
    if (rate.cat == 7) {
      rownames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)",
                                    "(0,R7)", "(1,R7)")
      colnames(rate.mat.index) <- c("(0,R1)", "(1,R1)",
                                    "(0,R2)", "(1,R2)", "(0,R3)", "(1,R3)", "(0,R4)",
                                    "(1,R4)", "(0,R5)", "(1,R5)", "(0,R6)", "(1,R6)",
                                    "(0,R7)", "(1,R7)")
    }
  }
  if (hrm == FALSE) {
    k = ntraits
    nl = 2
    if (ntraits == 1) {
      k <- 1
      nl <- nstates
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl, nl)
        tmp2 <- cbind(1:(nl^k), 1:(nl^k))
        index <- matrix(TRUE, nl^k, nl^k)
        diag(index) <- FALSE
        if (model == "ER") {
          np <- 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- nl * (nl - 1)/2
          sel <- col(rate.mat.index) < row(rate.mat.index)
          rate.mat.index <- t(rate.mat.index)
          rate.mat.index[sel] <- 1:np
          rate.mat.index[upper.tri(rate.mat.index)] = t(rate.mat.index)[upper.tri(rate.mat.index)]
        }
        if (model == "ARD") {
          np <- nl * (nl - 1)
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 2) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 0, 1, 1)
      vec.tmp2 <- c(0, 1, 0, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
      }
      matFINAL <- mat1 + mat2
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 4
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 4, 6)] <- rate.mat.index[index][c(3,
                                                                          5, 7, 8)] <- 1:np
        }
        if (model == "ARD") {
          np <- 8
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 3) {
      mat1 <- matrix(, nl^k, nl^k)
      mat2 <- matrix(, nl^k, nl^k)
      mat3 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 1, 0, 0, 1, 1, 0, 1)
      vec.tmp2 <- c(0, 0, 1, 0, 1, 0, 1, 1)
      vec.tmp3 <- c(0, 0, 0, 1, 0, 1, 1, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
        mat3[i, ] <- abs(vec.tmp3 - vec.tmp3[i])
      }
      matFINAL <- mat1 + mat2 + mat3
      if (is.character(model)) {
        rate.mat.index <- matrix(NA, nl^k, nl^k)
        if (model == "ER") {
          np <- 1
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
        if (model == "SYM") {
          np <- 12
          index <- matFINAL == 1
          rate.mat.index[index][c(1, 2, 3, 5, 6, 8, 9,
                                  11, 12, 15, 18, 21)] <- rate.mat.index[index][c(4,
                                                                                  7, 10, 13, 16, 14, 19, 17, 20, 22, 23, 24)] <- 1:np
        }
        if (model == "ARD") {
          np <- 24
          index <- matFINAL == 1
          rate.mat.index[index] <- 1:np
        }
      }
    }
    if (ntraits == 1) {
      rownames(rate.mat.index) <- as.character(1:nl)
      colnames(rate.mat.index) <- as.character(1:nl)
    }
    if (ntraits == 2) {
      rownames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)",
                                    "(1,1)")
      colnames(rate.mat.index) <- c("(0,0)", "(0,1)", "(1,0)",
                                    "(1,1)")
    }
    if (ntraits == 3) {
      rownames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)",
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)",
                                    "(1,1,1)")
      colnames(rate.mat.index) <- c("(0,0,0)", "(1,0,0)",
                                    "(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)",
                                    "(1,1,1)")
    }
  }
  return(rate.mat.index)
}
rate.par.drop <- function(rate.mat.index=NULL,drop.par=NULL){
  if(is.null(rate.mat.index)){
    stop("Rate matrix needed.  See mat.maker to create one.\n")
  }
  if(is.null(drop.par)){
    cat("No parameters indicated to drop.  Original matrix returned.\n")
    return(rate.mat.index)
  }
  if(max(rate.mat.index,na.rm=TRUE) < max(drop.par,na.rm=TRUE)){
    cat("Some parameters selected for dropping were not in the original matrix.\n")
  }
  drop.par <- unique(drop.par) # in case parameters listed more than once in drop vector
  drop.par <- drop.par[order(drop.par)]
  max <- max(rate.mat.index,na.rm=TRUE)
  for(drop.which in 1:length(drop.par)){
    drop.locs <- which(rate.mat.index == drop.par[drop.which],arr.ind=TRUE)
    rate.mat.index[drop.locs] <- NA
  }
  max <- max - length(drop.par)
  exclude <- which(is.na(rate.mat.index))
  rate.mat.index[-exclude] <- 1:max
  
  return(rate.mat.index)
}

#root function
root.obs <- function(prop, bins, cutoff =3.1){
  qq <- qnorm(prop, 0, 1)
  x <- pnorm(seq(-1*cutoff, cutoff, length.out=bins-1), qq,1)
  x <- c(x[2:(length(x)-1)],1) - c(0,x[2:(length(x)-1)])
  P <- c(0, x,0)
  P <- P/sum(P) #Make sure it sums to 1
  return(P)
}

#proportion function
prop <- function(data) {
  x = table(unlist(data[,2]))
  prop_1 <- x[1] / (x[1]+x[2])
  return(prop_1)
}

#BIC function
bic <- function(loglik, k, n) {
  ##where k is number of states 
  ##and n is number of taxa (up for debate lol)
  stopifnot(n > 0, k >= 0)
  bic_v <- -2*loglik + k*log(n)
  return(bic_v) 
} 

make_starting_values <- function(HMM_model,thresh_model,n_state) {
ratemat6 <- HMM_model$index.mat
ratemat6 <- ratemat6[c(rev(seq(1,2*n_state,2)),seq(2,2*n_state,2)),]
ratemat6  <-  ratemat6[,c(rev(seq(1,2*n_state,2)),seq(2,2*n_state,2))]
ratemat6[!is.na(ratemat6) & !is.na(thresh_model$solution)] <- thresh_model$solution[!is.na(ratemat6) & !is.na(thresh_model$solution)]
ratemat6[!is.na(ratemat6) & is.na(thresh_model$solution)] <- 1e-9
o <- order(c(rev(seq(1,2*n_state,2)),seq(2,2*n_state,2)))
ratemat6 <- ratemat6[o,]
ratemat6  <-  ratemat6[,o]
return(ratemat6)
}

#saving aic in dataframe 
put_df_aic <- function(data, fit, i){
  state <- ncol(fit$solution) 
  samp <- length(fit$phy$tip.label)
  fit_char <- deparse(substitute(fit)) 
  if(grepl('thresh', fit_char)) {
    mod_name <- "thresh"
  } else {
    mod_name <- "mk"
  }
  
  data[i,1] <- mod_name 
  data[i,2] <- state
  data[i,3] <- samp
  data[i,4] <- fit$AIC
  
  return(data[i,])
}

# basic Mk model 
## lists to contain the 200 sims for each Markov model and the threshold model
matrix_listARD <- replicate(n=8, expr=list())
matrix_listER <- replicate(n=8, expr=list())
matrix_listThreshER <- replicate(n=8, expr=list())
matrix_listThreshARD <- replicate(n=8, expr=list())
for (i in  1:8){
  #ARD
  m1 <- rate.mat.maker(ntraits = 1, nstates = i, model = "ARD", hrm = F)
  matrix_listARD[[i]]<- m1
  
  #ER
  m2 <- rate.mat.maker(ntraits = 1, nstates=i, model="ER", hrm = F)
  matrix_listER[[i]] <- m2
  
  #need to modify this for a threshold model
  m4 <- rate.mat.maker(ntraits = 1, nstates=i, model = "SYM", hrm=F)
  matrix_listThreshER[[i]] <- m4
}


mk_four <- matrix_listER[[4]]
mk_six <- matrix_listER[[6]]
mk_eight <- matrix_listER[[8]]

#mk_ten <- matrix(ncol = 10, nrow = 10,byrow = T, c(NA,1,2,3,4,5,6,7,8,9,
#                                                   10,NA,11,12,13,14,15,16,17,18,
#                                                  19,20,NA,21,22,23,24,25,26,27,
#                                                   28,29,30,NA,31,32,33,34,35,36,
#                                                   37,38,39,40,NA,42,43,44,45,46,
#                                                   47,48,49,50,51,NA,52,53,54,55,
#                                                   56,57,58,59,60,61,NA,62,63,64,
#                                                   65,66,67,68,69,70,71,NA,72,73,
#                                                   74,75,76,77,78,79,80,81,NA,82,
#                                                   83,84,85,86,87,88,89,90,91,NA))

mk_ten <- rate.mat.maker(ntraits =1, nstates =10, model = "ER", hrm=F)
mk_16 <- rate.mat.maker(ntraits=1, nstates = 16, model = "ER", hrm = F)
#mk_20 <- rate.mat.maker(ntraits=1, nstates = 20, model = "ER", hrm = F)

#hidden thresh models




two <- matrix_listThreshER[[2]]

## three


## four

four <- matrix(c(NA,1,NA,NA,NA,NA,1,NA,NA,1,NA,NA,NA,NA,1,NA), nrow=4, ncol=4)
rownames(four) <- c("1","2","3","4")
colnames(four) <- c("1","2","3","4")
matrix_listThreshER[[4]] <- four

## five -- don't really need this
five <- matrix(c(NA,1,NA,NA,NA,NA,NA,1,NA,NA,NA,1,NA,1,NA,NA,NA,1,NA,NA,NA,NA,NA,1,NA), nrow =5, ncol =5)
rownames(five) <- c("1","2","3","4","5")
colnames(five) <- c("1","2","3","4","5")
matrix_listThreshER[[5]] <- five

## six

six <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,NA,NA,NA,NA,NA,1,NA),
              nrow=6, ncol=6)
rownames(six) <- c("1","2","3","4","5","6")
colnames(six) <- c("1","2","3","4","5","6")
matrix_listThreshER[[6]] <- six

## seven
seven <- matrix(c(NA,1,NA,NA,NA,NA,NA,
                  1,NA,1,NA,NA,NA,NA
                  ,NA,1,NA,1,NA,NA,NA
                  ,NA,NA,1,NA,1,NA,NA,
                  NA,NA,NA,1,NA,1,NA,
                  NA,NA,NA,NA,1,NA,1,
                  NA,NA,NA,NA,NA,1,NA),
                ncol=7, nrow=7)
rownames(seven) <- c("1","2","3","4","5","6","7")
colnames(seven) <- c("1","2","3","4","5", "6","7")
matrix_listThreshER[[7]] <- seven


## eight
eight <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,
                  NA,NA,1,NA,NA,NA,NA, NA
                  ,NA,1,NA,1,NA,NA,NA,NA
                  ,NA,NA,1,NA,1,NA,NA,NA,
                  NA,NA,NA,1,NA,1,NA,NA,
                  NA,NA,NA,NA,1,NA,1,NA,
                  NA,NA,NA,NA,NA,1,NA,NA
                  ,NA,NA,NA,NA,NA,NA,1,NA
),
ncol=8, nrow=8)
rownames(eight) <- c("1","2","3","4","5","6","7","8")
colnames(eight) <- c("1","2","3","4","5", "6","7","8")
matrix_listThreshER[[8]] <- eight


# ER threshold matrices 

thresh_four <- matrix_listThreshER[[4]]

thresh_six <- matrix_listThreshER[[6]]

thresh_eight <- matrix_listThreshER[[8]]

thresh_ten <-  matrix(ncol = 10, nrow = 10,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

thresh_12 <- matrix(ncol = 12, nrow = 12,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))

thresh_14 <- matrix(ncol = 14, nrow = 14,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                      NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))




thresh_16 <-  matrix(ncol = 16, nrow = 16,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                       NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))



## SIMULATION 

sim_aic <- list()
everyting <- list()


# if errors out put in trycatch, 
# try with do first then dopar, then full

for(i in 1:nsim){
  all_aic <- matrix(nrow = 3, ncol = 7)
  count <- as.character(i)
  a0 <- runif(1,-1,1)
  
  #trees
  tree50 <- sim.bdtree(b=1,d = 0,n = 50)
  tree100 <- sim.bdtree(b=1, d=0, n=100)
  tree200 <- sim.bdtree(b=1,d=0, n=200)
  
  #data
  
  thresh50 <- fastBM(tree50, a = a0, sig2 =1)
  
  thresh100 <- fastBM(tree100, a=a0, sig2 = 1)

  thresh200 <- fastBM(tree200, a=a0, sig2 = 1)

  
  #splitting
  #50
  
  th50_df <- data.frame(thresh50)
  th100_df <- data.frame(thresh100)
  th200_df <- data.frame(thresh200)
  
  while(length((subset(th50_df, thresh50 < a0 )[,1])) == 0 || length((subset(th50_df, thresh50 > a0 )[,1])) == 0 
        || length((subset(th100_df, thresh100 > a0 )[,1])) == 0 || length((subset(th100_df, thresh100 < a0 )[,1])) == 0 ||
        length((subset(th200_df, thresh200 > a0 )[,1])) == 0 || length((subset(th200_df, thresh200 < a0 )[,1])) == 0 ) {
    
  a0 <- runif(1,-1,1)
  thresh50 <- fastBM(tree50, a = a0, sig2 =1)
  thresh100 <- fastBM(tree100, a=a0, sig2 = 1)
  thresh200 <- fastBM(tree200, a=a0, sig2 = 1)
  th50_df <- data.frame(thresh50)
  th100_df <- data.frame(thresh100)
  th200_df <- data.frame(thresh200)
  }
  
  
  
  th50_four <- case_when(
    th50_df[,1] > a0 ~ "1&2",
    th50_df[,1] < a0 ~ "3&4")
  # corr
  th50_four_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2")
  
  
  th50_six <- case_when(
    th50_df[,1] > a0 ~ "1&2&3",
    th50_df[,1] < a0 ~ "4&5&6"
  )
  # corr
  th50_six_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2"
  )
  
  th50_eight <- case_when(
    th50_df[,1] > a0 ~ "1&2&3&4",
    th50_df[,1] < a0 ~ "5&6&7&8"
  )
  #corr
  th50_eight_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2"
  )
  
  th50_ten <- case_when(
    th50_df[,1] > a0 ~ "1&2&3&4&5",
    th50_df[,1] < a0 ~ "6&7&8&9&10"
  )
  th50_ten_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2"
  )
  
  th50_12 <- case_when(
    th50_df[,1] > a0 ~ "1&2&3&4&5&6",
    th50_df[,1] < a0 ~ "7&8&9&10&11&12"
  )
  
  th50_12_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2"
  )
  
  th50_14 <- case_when(
    th50_df[,1] > a0 ~ "1&2&3&4&5&6&7",
    th50_df[,1] < a0 ~ "8&9&10&11&12&13&14"
  )
  
  th50_14_corr <- case_when(
    th50_df[,1] > a0 ~ "1",
    th50_df[,1] < a0 ~ "2"
  )
  
  
  th50_16 <- case_when(
    th50_df[,1] >a0 ~ "1&2&3&4&5&6&7&8",
    th50_df[,1] <a0 ~ "9&10&11&12&13&14&15&16"
  )
  
  th50_16_corr <- case_when(
    th50_df[,1] >a0 ~ "1",
    th50_df[,1] <a0 ~ "2"
  )
  
  th50_20 <- case_when(
    th50_df[,1] >a0 ~ "1&2&3&4&5&6&7&8&9&10",
    th50_df[,1] <a0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  th50_20_corr <- case_when(
    th50_df[,1] >a0 ~ "1",
    th50_df[,1] <a0 ~ "2"
  )
  
  
  
  th50_four <- data.frame(tree50$tip.label, th50_four)
  th50_four_corr <- data.frame(tree50$tip.label, th50_four_corr)
  th50_six <- data.frame(tree50$tip.label, th50_six)
  th50_six_corr <- data.frame(tree50$tip.label, th50_six_corr)
  th50_eight <- data.frame(tree50$tip.label, th50_eight)
  th50_eight_corr <- data.frame(tree50$tip.label, th50_eight_corr)
  th50_ten <- data.frame(tree50$tip.label, th50_ten)
  th50_ten_corr <- data.frame(tree50$tip.label, th50_ten_corr)
  th50_12 <- data.frame(tree50$tip.label, th50_12)
  th50_12_corr <- data.frame(tree50$tip.label, th50_12_corr)
  th50_14 <- data.frame(tree50$tip.label, th50_14)
  th50_14_corr <- data.frame(tree50$tip.label, th50_14_corr)
  th50_16 <- data.frame(tree50$tip.label, th50_16)
  th50_16_corr <- data.frame(tree50$tip.label, th50_16_corr)
  th50_20 <- data.frame(tree50$tip.label, th50_20)
  th50_20_corr <- data.frame(tree50$tip.label, th50_20_corr)
  
  #splitting
  #100

  th100_four <- case_when(
    th100_df[,1] > a0 ~ "1&2",
    th100_df[,1] < a0 ~ "3&4")
  # corr
  th100_four_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2")
  
  
  th100_six <- case_when(
    th100_df[,1] > a0 ~ "1&2&3",
    th100_df[,1] < a0 ~ "4&5&6"
  )
  # corr
  th100_six_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2"
  )
  
  th100_eight <- case_when(
    th100_df[,1] > a0 ~ "1&2&3&4",
    th100_df[,1] < a0 ~ "5&6&7&8"
  )
  #corr
  th100_eight_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2"
  )
  
  th100_ten <- case_when(
    th100_df[,1] > a0 ~ "1&2&3&4&5",
    th100_df[,1] < a0 ~ "6&7&8&9&10"
  )
  th100_ten_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2"
  )
  
  th100_12 <- case_when(
    th100_df[,1] > a0 ~ "1&2&3&4&5&6",
    th100_df[,1] < a0 ~ "7&8&9&10&11&12"
  )
  
  th100_12_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2"
  )
  
  th100_14 <- case_when(
    th100_df[,1] > a0 ~ "1&2&3&4&5&6&7",
    th100_df[,1] < a0 ~ "8&9&10&11&12&13&14"
  )
  
  th100_14_corr <- case_when(
    th100_df[,1] > a0 ~ "1",
    th100_df[,1] < a0 ~ "2"
  )
  
  th100_16 <- case_when(
    th100_df[,1] >a0 ~ "1&2&3&4&5&6&7&8",
    th100_df[,1] <a0 ~ "9&10&11&12&13&14&15&16"
  )
  
  th100_16_corr <- case_when(
    th100_df[,1] >a0 ~ "1",
    th100_df[,1] <a0 ~ "2"
  )
  
  th100_20 <- case_when(
    th100_df[,1] >a0 ~ "1&2&3&4&5&6&7&8&9&10",
    th100_df[,1] <a0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  th100_20_corr <- case_when(
    th100_df[,1] >a0 ~ "1",
    th100_df[,1] <a0 ~ "2"
  )
  
  
  th100_four <- data.frame(tree100$tip.label, th100_four)
  th100_four_corr <- data.frame(tree100$tip.label, th100_four_corr)
  th100_six <- data.frame(tree100$tip.label, th100_six)
  th100_six_corr <- data.frame(tree100$tip.label, th100_six_corr)
  th100_eight <- data.frame(tree100$tip.label, th100_eight)
  th100_eight_corr <- data.frame(tree100$tip.label, th100_eight_corr)
  th100_ten <- data.frame(tree100$tip.label, th100_ten)
  th100_ten_corr <- data.frame(tree100$tip.label, th100_ten_corr)
  th100_12 <- data.frame(tree100$tip.label, th100_12)
  th100_12_corr <- data.frame(tree100$tip.label, th100_12_corr)
  th100_14 <- data.frame(tree100$tip.label, th100_14)
  th100_14_corr <- data.frame(tree100$tip.label, th100_14_corr)
  th100_16 <- data.frame(tree100$tip.label, th100_16)
  th100_16_corr <- data.frame(tree100$tip.label, th100_16_corr)
  th100_20 <- data.frame(tree100$tip.label, th100_20)
  th100_20_corr <- data.frame(tree100$tip.label, th100_20_corr)
  
  #splitting
  #200
  
  th200_four <- case_when(
    th200_df[,1] > a0 ~ "1&2",
    th200_df[,1] < a0 ~ "3&4")
  # corr
  th200_four_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2")
  
  
  th200_six <- case_when(
    th200_df[,1] > a0 ~ "1&2&3",
    th200_df[,1] < a0 ~ "4&5&6"
  )
  # corr
  th200_six_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2"
  )
  
  th200_eight <- case_when(
    th200_df[,1] > a0 ~ "1&2&3&4",
    th200_df[,1] < a0 ~ "5&6&7&8"
  )
  #corr
  th200_eight_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2"
  )
  
  th200_ten <- case_when(
    th200_df[,1] > a0 ~ "1&2&3&4&5",
    th200_df[,1] < a0 ~ "6&7&8&9&10"
  )
  th200_ten_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2"
  )
  
  th200_12 <- case_when(
    th200_df[,1] > a0 ~ "1&2&3&4&5&6",
    th200_df[,1] < a0 ~ "7&8&9&10&11&12"
  )
  
  th200_12_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2"
  )
  
  th200_14 <- case_when(
    th200_df[,1] > a0 ~ "1&2&3&4&5&6&7",
    th200_df[,1] < a0 ~ "8&9&10&11&12&13&14"
  )
  
  th200_14_corr <- case_when(
    th200_df[,1] > a0 ~ "1",
    th200_df[,1] < a0 ~ "2"
  )
  
  th200_16 <- case_when(
    th200_df[,1] >a0 ~ "1&2&3&4&5&6&7&8",
    th200_df[,1] <a0 ~ "9&10&11&12&13&14&15&16"
  )
  
  th200_16_corr <- case_when(
    th200_df[,1] >a0 ~ "1",
    th200_df[,1] <a0 ~ "2"
  )
  
  th200_20 <- case_when(
    th200_df[,1] >a0 ~ "1&2&3&4&5&6&7&8&9&10",
    th200_df[,1] <a0 ~ "11&12&13&14&15&16&17&18&19&20"
  )
  th200_20_corr <- case_when(
    th200_df[,1] >a0 ~ "1",
    th200_df[,1] <a0 ~ "2"
  )
  
  
  th200_four <- data.frame(tree200$tip.label, th200_four)
  th200_four_corr <- data.frame(tree200$tip.label, th200_four_corr)
  th200_six <- data.frame(tree200$tip.label, th200_six)
  th200_six_corr <- data.frame(tree200$tip.label, th200_six_corr)
  th200_eight <- data.frame(tree200$tip.label, th200_eight)
  th200_eight_corr <- data.frame(tree200$tip.label, th200_eight_corr)
  th200_ten <- data.frame(tree200$tip.label, th200_ten)
  th200_ten_corr <- data.frame(tree200$tip.label, th200_ten_corr)
  th200_12 <- data.frame(tree200$tip.label, th200_12)
  th200_12_corr <- data.frame(tree200$tip.label, th200_12_corr)
  th200_14 <- data.frame(tree200$tip.label, th200_14)
  th200_14_corr <- data.frame(tree200$tip.label, th200_14_corr)
  th200_16 <- data.frame(tree200$tip.label, th200_16)
  th200_16_corr <- data.frame(tree200$tip.label, th200_16_corr)
  th200_20 <- data.frame(tree200$tip.label, th200_20)
  th200_20_corr <- data.frame(tree200$tip.label, th200_20_corr)
  
  
  # 4 thresh/Mk
  #start_four  <- matrix(c(0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0), nrow=4, ncol=4)
  ## 50
  #four50_Mk <- rayDISC(phy = tree50, data = th50_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th50_four), 4))
  
    # two hidden states 
  four50_corr <- corHMM(phy = tree50, data = th50_four_corr, rate.cat = 2, model = "ER", node.states = "none") 
  four50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th50_four), 4), node.states = "none")
  
  n_state <- 2
  
  # makes starting values at least related to threshold model, to give corhmm a fighting chance 
  
  ratemat4 <- make_starting_values(four50_corr, four50_thresh, n_state)
  p4 <- ratemat4[!is.na(as.vector(ratemat4))]
  ratemat4[!is.na(ratemat4)] <- 1:8
  ratemat_final4 <- ratemat4
  ratemat_final4[!is.na(ratemat4)] <- p4
  p4[p4 <= 1e-9] <- 1e-9
  
    #p4 <- as.vector(four50_corr$solution)
    #p4 <- na.omit(p4)
    #ratemat4 <- four50_corr$index.mat
    #ratemat4[!is.na(ratemat4)] <- 1:8
    #ratemat_final4 <- ratemat4
    #ratemat_final4[!is.na(ratemat4)] <- p4
  
 four50_corr <- corHMM(phy = tree50, data = th50_four_corr, rate.cat = 2, rate.mat = ratemat4, ip = p4, node.states = "none") 
    
  
  ## 200
  #four100_Mk <- rayDISC(phy= tree100, data = th100_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th100_four), 4))
    
  four100_corr <- corHMM(phy = tree100, data = th100_four_corr, rate.cat  = 2, model = "ER", node.states = "none")
  four100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th100_four), 4), node.states = "none")
  
    ratemat4 <- make_starting_values(four100_corr, four100_thresh, n_state)
    p4 <- ratemat4[!is.na(as.vector(ratemat4))]
    ratemat4[!is.na(ratemat4)] <- 1:8
    ratemat_final4 <- ratemat4
    ratemat_final4[!is.na(ratemat4)] <- p4
    p4[p4 <= 1e-9] <- 1e-9

  four100_corr <- corHMM(phy = tree100, data = th100_four_corr, rate.cat = 2, rate.mat = ratemat4, ip = p4, node.states = "none") 

  ## 200
  #four200_Mk <- rayDISC(phy = tree200, data = th200_four, rate.mat = mk_four, root.p = root.obs(prop = prop(th200_four), 4))
      
  four200_corr <- corHMM(phy = tree200, data = th200_four_corr, rate.cat = 2, model = "ER", node.states = "none")
  four200_thresh <- rayDISC(model = "ER", phy = tree200, data = th200_four, rate.mat = thresh_four, root.p = root.obs(prop = prop(th200_four), 4), node.states = "none")
  
    ratemat4 <- make_starting_values(four200_corr, four200_thresh, n_state)
    p4 <- ratemat4[!is.na(as.vector(ratemat4))]
    ratemat4[!is.na(ratemat4)] <- 1:8
    ratemat_final4 <- ratemat4
    ratemat_final4[!is.na(ratemat4)] <- p4
    p4[p4 <= 1e-9] <- 1e-9
    
  four200_corr <- corHMM(phy = tree200, data = th200_four_corr, rate.cat = 2, rate.mat = ratemat4, ip = p4, node.states = "none")

  # 6 thresh/Mk
  
  ## 50
  #six50_Mk <- rayDISC(phy = tree50, data = th50_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th50_six), 6))
  n_state <- 3
  six50_corr <- corHMM(phy = tree50, data = th50_six_corr, rate.cat = 3, model = "ER", node.states = "none")
  six50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_six, rate.mat = thresh_six, root.p =  root.obs(prop = prop(th50_six), 6), node.states = "none")
  
    ratemat6 <- make_starting_values(six50_corr, six50_thresh, n_state)
    p6 <- ratemat6[!is.na(as.vector(ratemat6))]
    ratemat6[!is.na(ratemat6)] <- 1:18
    ratemat_final6 <- ratemat6
    ratemat_final6[!is.na(ratemat6)] <- p6
    p6[p6 <= 1e-9] <- 1e-9
  
  six50_corr <- corHMM(phy = tree50, data = th50_six_corr, rate.cat = 3, rate.mat = ratemat6, ip = p6, node.states =  "none") 

  ## 100
  #six100_Mk <- rayDISC(phy= tree100, data = th100_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th100_six), 6))
  n_state <- 3
  six100_corr <- corHMM(phy = tree100, data = th100_six_corr, rate.cat = 3, model = "ER", node.states = "none")
  six100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_six, rate.mat = thresh_six, root.p = root.obs(prop = prop(th100_six), 6), node.states = "none")

    ratemat6 <- make_starting_values(six100_corr, six100_thresh, n_state)
    p6 <- ratemat6[!is.na(as.vector(ratemat6))]
    ratemat6[!is.na(ratemat6)] <- 1:18
    ratemat_final6 <- ratemat6
    ratemat_final6[!is.na(ratemat6)] <- p6
    p6[p6 <= 1e-9] <- 1e-9

  six100_corr <- corHMM(phy = tree100, data = th100_six_corr, rate.cat = 3, rate.mat = ratemat6, ip = p6, node.states = "none") 
  
  
  ## 200
  #six200_Mk <- rayDISC(phy = tree200, data = th200_six, rate.mat = mk_six, root.p = root.obs(prop = prop(th200_six), 6))
  
  six200_corr <- corHMM(phy = tree200, data = th200_six_corr, rate.cat = 3, model = "ER", node.states = "none")
  
  six200_thresh <- rayDISC(model = "ER", phy = tree200, data = th200_six, rate.mat = thresh_six, root.p = root.obs(prop = prop(th200_six), 6), node.states = "none")
  
    ratemat6 <- make_starting_values(six200_corr, six200_thresh, 3)
    p6 <- ratemat6[!is.na(as.vector(ratemat6))]
    ratemat6[!is.na(ratemat6)] <- 1:18
    ratemat_final6 <- ratemat6
    ratemat_final6[!is.na(ratemat6)] <- p6
    p6[p6 <= 1e-9] <- 1e-9

  six200_corr <- corHMM(phy = tree200, data = th200_six_corr, rate.cat = 3, rate.mat = ratemat6, ip = p6, node.states = "none") 
  
  # 8 thresh/mk/
  
  ## 50
  #eight50_Mk <- rayDISC(phy = tree50, data = th50_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th50_eight), 8))
  
  eight50_corr <- corHMM(phy = tree50, data = th50_eight_corr, rate.cat = 4, model = "ER", node.states = "none")
  eight50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th50_eight), 8), node.states = "none")
   
     n_state <- 4
  
    ratemat8 <- make_starting_values(eight50_corr, eight50_thresh, n_state)
    p8 <- ratemat8[!is.na(as.vector(ratemat8))]
    ratemat8[!is.na(ratemat8)] <- 1:16
    ratemat_final8 <- ratemat8
    ratemat_final8[!is.na(ratemat8)] <- p8
    p8[p8 <= 1e-9] <- 1e-9
  
      #p8 <- as.vector(eight50_corr$solution)
      #p8 <- na.omit(p8)
      #ratemat8 <- eight50_corr$index.mat
      #ratemat8[!is.na(ratemat8)] <- 1:16
      #ratemat_final8 <- ratemat8
      #ratemat_final8[!is.na(ratemat8)] <- p8
  
  eight50_corr <- corHMM(phy = tree50, data = th50_eight_corr, rate.cat = 4, rate.mat = ratemat8, ip = p8, node.states = "none") 
  
  ## 100
  #eight100_Mk <- rayDISC(phy= tree100, data = th100_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th100_eight), 8))
  
  eight100_corr <- corHMM(phy = tree100, data = th100_eight_corr, rate.cat = 4, model = "ER", node.states = "none")
  eight100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th100_eight), 8), node.states = "none")
    
    ratemat8 <- make_starting_values(eight100_corr, eight100_thresh, n_state)
    p8 <- ratemat8[!is.na(as.vector(ratemat8))]
    ratemat8[!is.na(ratemat8)] <- 1:16
    ratemat_final8 <- ratemat8
    ratemat_final8[!is.na(ratemat8)] <- p8
    p8[p8 <= 1e-9] <- 1e-9
  
  eight100_corr <- corHMM(phy = tree100, data = th100_eight_corr, rate.cat = 4, rate.mat = ratemat8, ip = p8, node.states = "none") 
    
  ## 200
  #eight200_Mk <- rayDISC(phy = tree200, data = th200_eight, rate.mat = mk_eight, root.p = root.obs(prop = prop(th200_eight), 8))
  
  eight200_corr <- corHMM(phy = tree200, data = th200_eight_corr, rate.cat = 4, model = "ER", node.states = "none")
    
  eight200_thresh <- rayDISC(model = "ER", phy = tree200, data = th200_eight, rate.mat = thresh_eight, root.p = root.obs(prop = prop(th200_eight), 8), node.states = "none")
    
    ratemat8 <- make_starting_values(eight200_corr, eight200_thresh, n_state)
    p8 <- ratemat8[!is.na(as.vector(ratemat8))]
    ratemat8[!is.na(ratemat8)] <- 1:16
    ratemat_final8 <- ratemat8
    ratemat_final8[!is.na(ratemat8)] <- p8
    p8[p8 <= 1e-9] <- 1e-9
    
    #p8 <- as.vector(eight200_corr$solution)
    #p8 <- na.omit(p8)
    #ratemat8 <- eight200_corr$index.mat
    #ratemat8[!is.na(ratemat8)] <- 1:16
    #ratemat_final8 <- ratemat8
    #ratemat_final8[!is.na(ratemat8)] <- p8
    
  eight200_corr <- corHMM(phy = tree200, data = th200_four_corr, rate.cat = 4, rate.mat = ratemat8, ip = p8, node.states = "none") 

  # 10 thresh/mk
  
  ## 50
  #ten50_Mk <- rayDISC(phy = tree50, data = th50_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th50_ten), 10))
  n_state <- 5
  ten50_corr <- corHMM(phy = tree50, data = th50_ten_corr, rate.cat = 5, model = "ER", node.states = "none")
  ten50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th50_ten), 10), node.states = "none")
    
    #p10 <- as.vector(ten50_corr$solution)
    #p10 <- na.omit(p10)
    #ratemat10 <- ten50_corr$index.mat
    # this needs to repeat
    #ratemat10[!is.na(ratemat10)] <-c(1:20, 1:20, 1:10)
    #ratemat_final10 <- ratemat10
    #ratemat_final10[!is.na(ratemat10)] <- p10
    
      ratemat10 <- make_starting_values(ten50_corr, ten50_thresh, n_state)
      p10 <- ratemat10[!is.na(as.vector(ratemat10))]
      ratemat10[!is.na(ratemat10)] <- 1:50
      ratemat_final10 <- ratemat10
      ratemat_final10[!is.na(ratemat10)] <- p10
      p10[p10 <= 1e-9] <- 1e-9
    
  ten50_corr <- corHMM(phy = tree50, data = th50_ten_corr, rate.cat = 5, rate.mat = ratemat10, ip = p10, node.states = "none") 
    
  ## 100
  #ten100_Mk <- rayDISC(phy= tree100, data = th100_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th100_ten), 10))
  
  
  ten100_corr <- corHMM(phy = tree100, data = th100_ten_corr, rate.cat = 5, model = "ER", node.states = "none")
    
  ten100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th100_ten), 10), node.states = "none")
  
    #p10 <- as.vector(ten100_corr$solution)
    #p10 <- na.omit(p10)
    #ratemat10 <- ten100_corr$index.mat
    #ratemat10[!is.na(ratemat10)] <- c(1:20, 1:20, 1:10)
    #ratemat_final10 <- ratemat10
    #ratemat_final10[!is.na(ratemat10)] <- p10
    
    ratemat10 <- make_starting_values(ten100_corr, ten100_thresh, n_state)
    p10 <- ratemat10[!is.na(as.vector(ratemat10))]
    ratemat10[!is.na(ratemat10)] <- 1:50
    ratemat_final10 <- ratemat10
    ratemat_final10[!is.na(ratemat10)] <- p10
    p10[p10 <= 1e-9] <- 1e-9
    
  ten100_corr <- corHMM(phy = tree100, data = th100_ten_corr, rate.cat = 5, rate.mat = ratemat10, ip = p10, node.states = "none") 
  
  ## 200
  #ten200_Mk <- rayDISC(phy = tree200, data = th200_ten, rate.mat = mk_ten, root.p = root.obs(prop = prop(th200_ten), 10))
  
  ten200_corr <- corHMM(phy = tree200, data = th200_ten_corr, rate.cat = 5, model = "ER", node.states = "none")
  ten200_thresh <- rayDISC(model = "ER", phy = tree200, data = th200_ten, rate.mat = thresh_ten, root.p = root.obs(prop = prop(th200_ten), 10), node.states = "none")
  
    ratemat10 <- make_starting_values(ten200_corr, ten200_thresh, n_state)
    p10 <- ratemat10[!is.na(as.vector(ratemat10))]
    ratemat10[!is.na(ratemat10)] <- 1:50
    ratemat_final10 <- ratemat10
    ratemat_final10[!is.na(ratemat10)] <- p10
    p10[p10 <= 1e-9] <- 1e-9
  
  ten200_corr <- corHMM(phy = tree200, data = th200_ten_corr, rate.cat = 5, rate.mat = ratemat10, ip = p10, node.states ="none") 
  
  
  
  # 12 thresh/corhmm
  ## 50 
  twelve50_corr <- corHMM(phy = tree50, data = th50_12_corr, rate.cat = 6, model = "ER", node.states = "none")
  
  twelve50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_12, rate.mat = thresh_12, root.p = root.obs(prop = prop(th50_12), 12), node.states = "none")
  n_state <- 6
  ratemat12 <- make_starting_values(twelve50_corr, twelve50_thresh, n_state)
  p12 <- ratemat12[!is.na(as.vector(ratemat12))]
  ratemat12[!is.na(ratemat12)] <- 1:n_state^2
  ratemat_final12 <- ratemat12
  ratemat_final12[!is.na(ratemat12)] <- p12
  p12[p12 <= 1e-9] <- 1e-9
  
  twelve50_corr <- corHMM(phy = tree50, data = th50_12_corr, rate.cat = 6, rate.mat = ratemat12, ip = p12, node.states = "none") 
  
  ## 100
  
  twelve100_corr <- corHMM(phy = tree100, data = th100_12_corr, rate.cat = 6, model = "ER", node.states = "none")
  
  twelve100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_12, rate.mat = thresh_12, root.p = root.obs(prop = prop(th100_12), 12), node.states = "none")
  n_state <- 6
    ratemat12 <- make_starting_values(twelve100_corr, twelve100_thresh, n_state)
    p12 <- ratemat12[!is.na(as.vector(ratemat12))]
    ratemat12[!is.na(ratemat12)] <- 1:n_state^2
    ratemat_final12 <- ratemat12
    ratemat_final12[!is.na(ratemat12)] <- p12
    p12[p12 <= 1e-9] <- 1e-9
  
  twelve100_corr <- corHMM(phy = tree100, data = th100_12_corr, rate.cat = 6, rate.mat = ratemat12, ip = p12, node.states = "none") 
  
  
  ## 200 
  
  twelve200_corr <- corHMM(phy = tree200, data = th200_12_corr, rate.cat = 6, model = "ER", node.states = "none")
  
  twelve200_thresh <- rayDISC(model = "ER", phy= tree200, data = th200_12, rate.mat = thresh_12, root.p = root.obs(prop = prop(th200_12), 12), node.states = "none")
  n_state <- 6
  ratemat12 <- make_starting_values(twelve200_corr, twelve200_thresh, n_state)
  p12 <- ratemat12[!is.na(as.vector(ratemat12))]
  ratemat12[!is.na(ratemat12)] <- 1:n_state^2
  ratemat_final12 <- ratemat12
  ratemat_final12[!is.na(ratemat12)] <- p12
  p12[p12 <= 1e-9] <- 1e-9
  
  twelve200_corr <- corHMM(phy = tree200, data = th200_12_corr, rate.cat = 6, rate.mat = ratemat12, ip = p12, node.states = "none") 
  
  

  # 14 thresh/corhmm
  
  ## 50 
  fourteen50_corr <- corHMM(phy = tree50, data = th50_14_corr, rate.cat = 7, model = "ER", node.states = "none")
  
  fourteen50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th50_14), 14), node.states = "none")
  n_state <- 7
  ratemat14 <- make_starting_values(fourteen50_corr, fourteen50_thresh, n_state)
  p14 <- ratemat14[!is.na(as.vector(ratemat14))]
  ratemat14[!is.na(ratemat14)] <- 1:n_state^2
  ratemat_final14 <- ratemat14
  ratemat_final14[!is.na(ratemat14)] <- p14
  p14[p14 <= 1e-9] <- 1e-9
  
  fourteen50_corr <- corHMM(phy = tree50, data = th50_14_corr, rate.cat = 7, rate.mat = ratemat14, ip = p14, node.states = "none") 
  
  ## 100
  
  fourteen100_corr <- corHMM(phy = tree100, data = th100_14_corr, rate.cat = 7, model = "ER", node.states = "none")
  
  fourteen100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th100_14), 14), node.states = "none")
  n_state <- 7
  ratemat14 <- make_starting_values(fourteen100_corr, fourteen100_thresh, n_state)
  p14 <- ratemat14[!is.na(as.vector(ratemat14))]
  ratemat14[!is.na(ratemat14)] <- 1:n_state^2
  ratemat_final14 <- ratemat14
  ratemat_final14[!is.na(ratemat14)] <- p14
  p14[p14 <= 1e-9] <- 1e-9
  
  fourteen100_corr <- corHMM(phy = tree100, data = th100_14_corr, rate.cat = 7, rate.mat = ratemat14, ip = p14, node.states = "none") 
  
  
  ## 200 
  
  fourteen200_corr <- corHMM(phy = tree200, data = th200_14_corr, rate.cat = 7, model = "ER", node.states = "none")
  
  fourteen200_thresh <- rayDISC(model = "ER", phy= tree200, data = th200_14, rate.mat = thresh_14, root.p = root.obs(prop = prop(th200_14), 14), node.states = "none")
  n_state <- 7
  ratemat14 <- make_starting_values(fourteen200_corr, fourteen200_thresh, n_state)
  p14 <- ratemat14[!is.na(as.vector(ratemat14))]
  ratemat14[!is.na(ratemat14)] <- 1:n_state^2
  ratemat_final14 <- ratemat14
  ratemat_final14[!is.na(ratemat14)] <- p14
  p14[p14 <= 1e-9] <- 1e-9
  
  fourteen200_corr <- corHMM(phy = tree200, data = th200_14_corr, rate.cat = 7, rate.mat = ratemat14, ip = p14, node.states = "none") 
  
  
  # 16 thresh/corhmm
  
  ## 50
  
  sixteen50_corr <- corHMM(phy = tree50, data = th50_16_corr, rate.cat = 8, model = "ER", node.states = "none")
  sixteen50_thresh <- rayDISC(model = "ER", phy= tree50, data = th50_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th50_16), 16), node.states = "none")
  
  n_state <- 8
  
    ratemat16 <- make_starting_values(sixteen50_corr, sixteen50_thresh, n_state)
    p16 <- ratemat16[!is.na(as.vector(ratemat16))]
    ratemat16[!is.na(ratemat16)] <- 1:128
    ratemat_final16 <- ratemat16
    ratemat_final16[!is.na(ratemat16)] <- p16
    p16[p16 <= 1e-9] <- 1e-9
    
  
  sixteen50_corr <- corHMM(phy = tree50, data = th50_16_corr, rate.cat = 8, rate.mat = ratemat16, ip = p16, node.states = "none") 
  
  
  ## 100

   sixteen100_corr <- corHMM(phy = tree100, data = th100_16_corr, rate.cat = 8, model = "ER", node.states = "none")
   sixteen100_thresh <- rayDISC(model = "ER", phy= tree100, data = th100_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th100_16), 16), node.states = "none")
     n_state <- 8
     ratemat16 <- make_starting_values(sixteen100_corr, sixteen100_thresh, n_state)
     p16 <- ratemat16[!is.na(as.vector(ratemat16))]
     ratemat16[!is.na(ratemat16)] <- 1:128
     ratemat_final16 <- ratemat16
     ratemat_final16[!is.na(ratemat16)] <- p16
     p16[p16 <= 1e-9] <- 1e-9
   
   sixteen100_corr <- corHMM(phy = tree100, data = th100_16_corr, rate.cat = 8, rate.mat = ratemat16, ip = p16, node.states = "none") 
   

  ## 200
    sixteen200_corr <- corHMM(phy = tree200, data = th200_16_corr, rate.cat = 8, model = "ER", node.states = "none")
    sixteen200_thresh <- rayDISC(model = "ER", phy = tree200, data = th200_16, rate.mat = thresh_16, root.p = root.obs(prop = prop(th200_16),16), node.states = "none")
    
    ratemat16 <- make_starting_values(sixteen200_corr, sixteen200_thresh, n_state)
      p16 <- ratemat16[!is.na(as.vector(ratemat16))]
      ratemat16[!is.na(ratemat16)] <- 1:128
      ratemat_final16 <- ratemat16
      ratemat_final16[!is.na(ratemat16)] <- p16
      p16[p16 <= 1e-9] <- 1e-9
  
    
    sixteen200_corr <- corHMM(phy = tree200, data = th200_16_corr, rate.cat = 8, rate.mat = ratemat16, ip = p16, node.states = "none") 

  #delta AIC and BIC
  
  #four AIC
  all_aic[1,1] <- delta50_four <- four50_thresh$AIC - four50_corr$AIC
  all_aic[2,1] <- delta100_four <- four100_thresh$AIC - four100_corr$AIC
  all_aic[3,1] <- delta200_four <- four200_thresh$AIC - four200_corr$AIC

  #six AIC
  all_aic[1,2] <- delta50_six <- six50_thresh$AIC - six50_corr$AIC
  all_aic[2,2] <- delta100_six <- six100_thresh$AIC - six100_corr$AIC
  all_aic[3,2] <- delta200_six <- six200_thresh$AIC - six200_corr$AIC

  #eight AIC
  all_aic[1,3] <- delta50_eight <- eight50_thresh$AIC - eight50_corr$AIC
  all_aic[2,3] <- delta100_eight <- eight100_thresh$AIC - eight100_corr$AIC
  all_aic[3,3] <- delta200_eight <- eight200_thresh$AIC - eight200_corr$AIC
  

  #ten AIC
  all_aic[1,4] <- delta50_ten <- ten50_thresh$AIC - ten50_corr$AIC
  all_aic[2,4] <- delta100_ten <- ten100_thresh$AIC - ten100_corr$AIC
  all_aic[3,4] <- delta200_ten <- ten200_thresh$AIC - ten200_corr$AIC
  
  #twelve AIC
  all_aic[1,5] <- delta50_twelve <- twelve50_thresh$AIC - twelve50_corr$AIC
  all_aic[2,5] <- delta100_twelve <- twelve100_thresh$AIC - twelve100_corr$AIC
  all_aic[3,5] <- delta200_twelve <- twelve200_thresh$AIC - twelve200_corr$AIC
  
  #fourteen AIC
  all_aic[1,6] <- delta50_fourteen <- fourteen50_thresh$AIC - fourteen50_corr$AIC
  all_aic[2,6] <- delta100_fourteen <- fourteen100_thresh$AIC - fourteen100_corr$AIC
  all_aic[3,6] <- delta200_fourteen <- fourteen200_thresh$AIC - fourteen200_corr$AIC
  
  #sixteen AIC
  all_aic[1,7] <- delta50_sixteen <- sixteen50_thresh$AIC - sixteen50_corr$AIC
  all_aic[2,7] <- delta100_sixteen <- sixteen100_thresh$AIC - sixteen100_corr$AIC
  all_aic[3,7] <- delta200_sixteen <- sixteen200_thresh$AIC - sixteen200_corr$AIC
  
  # SAVING EACH SIM
  sim_aic[[i]] <- all_aic
}


## AIC 
all_sims <- sim_aic
setwd("~/heatmap_results/")
saveRDS(all_sims, file = "heatmap3_aic_results.rds")
