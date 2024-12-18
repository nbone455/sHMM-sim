#Empirical example
# corhmm 4 state vs thresh 4 state 
# whats the threshold trait here?
library(ape)
# Lets assume you have a trait (presence or absence of a third molar) that you 
# assume is a threshold trait -- an actually continuous trait that relies on 
setwd("../data/")
#prim_tree <- read.nexus("tree.nex")
prim_tree <- read.nexus("../data/prim_10k.nex")

# primate dataset
#dat <- read.csv("Plavcan_Lower Molar_All.csv")
# try with these, use the ratio as cutoff 
dat <- read.csv("Plavcan_Lower Molar_All.csv")
head(dat)
new_dat <- dat
new_dat[is.na(new_dat)] <- 0

# modifying species names 


#lets try this loop still 
  pres <- list()
  names <- list()
for(i in 1:length(new_dat$Row.names)) {
  if(is.na(new_dat$m3.md[i]) == TRUE || new_dat$m3.md[i] == 0){
    pres[[i]] <- 0
    names[[i]] <- new_dat$Row.names[i]
  } else {
    pres[[i]] <- 1 
    names[[i]] <- new_dat$Row.names[i]
    
  }
  
}
p_dat <- data.frame(unlist(pres))
p_dat$species <- unlist(names)
colnames(p_dat) <- c("pres.m3", "taxa")

# need to ratio the same spp. // dont need this

#count_ratio <- list()
#sp <- unique(p_dat$taxa)
#dplyr way
#library(dplyr)
#for(j in 1:length(sp)) {
#  t <- dplyr::filter(p_dat, taxa == sp[j]) 
#      if(0 %in% t$pres.m3 == TRUE){
#        count_ratio[j] <- length(which(t$pres.m3 == 0)) / length(t$pres.m3)
#        print(j)
#        } else {
#          count_ratio[j] <- 1
      #print("one")
#      }
#}
#counts <- data.frame(unlist(count_ratio))
#counts$taxa <- unlist(sp)
#colnames(counts) <- c("ratio", "taxa")

# discrete version

#d_dat <- case_when(
#  counts$ratio < 1 ~ 0,
#  counts$ratio == 1 ~ 1 
#)
#d_dat <- data.frame(counts$taxa, d_dat)
#colnames(d_dat) <- c("taxa","pres")

# three hidden states minimum 
thresh_dat <- case_when(
  p_dat$pres.m3 == 1 ~ "4&5&6",
  p_dat$pres.m3 == 0 ~ "1&2&3"
)
thresh_dat <- data.frame(p_dat$taxa, thresh_dat)

# four hidden states 
thresh_dat4 <- case_when(
  p_dat$pres.m3 == 1 ~ "5&6&7&8",
  p_dat$pres.m3 == 0 ~ "1&2&3&4"
)
thresh_dat4 <- data.frame(p_dat$taxa, thresh_dat4)

thresh_dat5 <- case_when(
  p_dat$pres.m3 == 1 ~ "6&7&8&9&10",
  p_dat$pres.m3 == 0 ~ "1&2&3&4&5"
)
thresh_dat5 <- data.frame(p_dat$taxa, thresh_dat5)


# removing function

#sp <- sort(sp)
#drop_subspecies <- function(x){
#y <- strsplit(x, "_")[[1]]
#z <- paste(y[1], y[2], sep="_")
#return(z)
#}

#new <- list()
# update this at a later time 
#some_trim <- c(5,7,8,9,25,11,15,25,26,27,30)
#sp_drop <- sp[some_trim]

#    for (k in 1:length(some_trim)) {
#      new[k] <- drop_subspecies(sp_drop[k])
#    }


#some_change <- c(64, )
# trying to replace subspecies with general sp.
#place <- which(d_dat$taxa %in% sp_drop) 
#new <- unlist(new)
#new <- new[-c(4,8)]
#d_dat$taxa[place] <- new


# matching tree & analysis 

# flipping columns 
p_dat <- p_dat[,c(2,1)]

library(treeplyr)
td <- make.treedata(prim_tree, p_dat)
#td_thresh <- make.treedata(prim_tree, thresh_dat)
corr_2 <- corHMM::corHMM(phy = td$phy, data = p_dat, rate.cat = 2)
corr_2ER <- corHMM::corHMM(phy = td$phy, data = p_dat, rate.cat = 2, model = "ER")

corr_1 <- corHMM(phy = td$phy, data = p_dat, rate.cat = 1)


corr_8 <- corHMM(phy = td$phy, data = p_dat, rate.cat = 4)

corr_10_ER <- corHMM(phy = td$phy, data = p_dat, rate.cat = 5, model = "ER")
corr_10_base <- corHMM(phy = td$phy, data = p_dat, rate.cat = 5)
## thresh model
three <- matrix( c(NA,1,NA,1,NA,1,NA,1,NA), nrow = 3, ncol= 3)
rownames(three) <- c("1","2","3")
colnames(three) <- c("1","2","3")
four <- matrix(c(NA,1,NA,NA,NA,NA,1,NA,NA,1,NA,NA,NA,NA,1,NA), nrow=4, ncol=4)
rownames(four) <- c("1","2","3","4")
colnames(four) <- c("1","2","3","4")
six <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,1,NA,NA,NA,NA,1,NA,NA,NA,NA,NA,NA,1,NA),
              nrow=6, ncol=6)
rownames(six) <- c("1","2","3","4","5","6")
colnames(six) <- c("1","2","3","4","5","6")

eight <- matrix(c(NA,1,NA,NA,NA,NA,NA,NA,
                       NA,NA,3,NA,NA,NA,NA, NA
                       ,NA,4,NA,5,NA,NA,NA,NA
                       ,NA,NA,6,NA,7,NA,NA,NA,
                       NA,NA,NA,8,NA,9,NA,NA,
                       NA,NA,NA,NA,10,NA,11,NA,
                       NA,NA,NA,NA,NA,12,NA,NA,
                       NA,NA,NA,NA,NA,NA,14,NA
),
ncol=8, nrow=8)
rownames(eight) <- c("1","2","3","4","5","6","7","8")
colnames(eight) <- c("1","2","3","4","5", "6","7","8")

ten <-  matrix(ncol = 10, nrow = 10,byrow = T, c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                                        1,NA,1,NA,NA,NA,NA,NA,NA,NA,
                                                        NA,1,NA,1,NA,NA,NA,NA,NA,NA,
                                                        NA,NA,1,NA,1,NA,NA,NA,NA,NA,
                                                        NA,NA,NA,1,NA,1,NA,NA,NA,NA,
                                                        NA,NA,NA,NA,1,NA,1,NA,NA,NA,
                                                        NA,NA,NA,NA,NA,1,NA,1,NA,NA,
                                                        NA,NA,NA,NA,NA,NA,1,NA,1,NA,
                                                        NA,NA,NA,NA,NA,NA,NA,1,NA,1,
                                                        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))


#corr_thresh <- corHMM(phy = td$phy, d_dat, rate.mat = four,rate.cat = 2)
 

# I Guess I need to fix the roots?
root.obs <- function(prop, bins, cutoff=3.1){
  qq <- qnorm(prop, 0, 1)
  x <- pnorm(seq(-1*cutoff, cutoff, length.out=bins-1), qq,1)
  x <- c(x[2:(length(x)-1)],1) - c(0,x[2:(length(x)-1)])
  P <- c(0, x,0)
  P <- P/sum(P) #Make sure it sums to 1
  return(P)
}
#proportion fx
prop <- function(data) {
  x = table(unlist(data[,2]))
  prop_1 <- x[1] / (x[1]+x[2])
  return(prop_1)
}

thresh_8 <- rayDISC(phy = td$phy, thresh_dat4, rate.mat = eight, root.p = root.obs(prop = prop(thresh_dat4), 8))

thresh_10 <- rayDISC(phy = td$phy, thresh_dat5, rate.mat = ten, root.p = root.obs(prop= prop(thresh_dat5), 10), model = "ER")




ray_ARD <- rayDISC(phy=td$phy, thresh_dat, rate)



#ray_thresh$phy$edge.length <- unlist(ray_thresh$phy$edge.length)
#paintBranches(ray_thresh$phy, ray_thresh$phy$edge.length, state = ray_thresh$states)

plotRECON(ray_thresh$phy, ray_thresh$states, piecolors = 
c("yellow", "cornflowerblue", "orange", "dark blue"), cex = 0.40, pie.cex = 0.25, height = 50, width = 12, 
show.tip.label = TRUE, label.offset = 0.4, title = "threshold molar",file = "threshEX.pdf", type = "fan")

#pure discrete plot
#edgey <- unlist(corr_two$phy$edge.length)
#paintBranches(corr_two$phy, edge = corr_two$phy$edge.length, state = corr_two$states)

corr_two_plot <- plotRECON(corr_two$phy, corr_two$states, piecolors = 
            c( "cornflowerblue", "orange"), cex = 0.40, pie.cex = 0.25, height = 50, width = 12, 
          show.tip.label = TRUE, label.offset = 0.4, title = "discrete molar",file = "discreteEX.pdf", type = "fan")

corr_four_plot <- plotRECON(corr_d$phy, corr_d$states, piecolors = 
                             c("yellow", "cornflowerblue", "orange", "dark blue"), cex = 0.40, pie.cex = 0.25, height = 50, width = 12, 
                           show.tip.label = TRUE, label.offset = 0.4, title = "corr 4-state molar",file = "corrEX.pdf", type = "fan")



# hidden state liability vs ratio (m2 / m1) vs (m3 / m1)

# liabilities are the posterior probabilities of being in one of the hidden states 


#subestting the orignal dataset 
unqNames <- unique(td$phy$tip.label)

just <- subset(new_dat, Row.names %in% unqNames)

plot(just$ratio_m3, thresh_10$states)

liabs_old <- thresh_10$tip.states %*% 1:10
liabs <- thresh_10$states[,(c(1,3,4,5,6,7,8,9,2))] %*% 1:9
l_tips <- liabs[1:length(td$phy$tip.label)]
liabs <- thresh_10$tip.states %*% 1:10

plot(thresh_10$states[,(c(1,3,4,5,6,7,8,9,2))] %*% 1:9, just$ratio_m3)

# this is for presence or absence

# mcmc threshold 
library(MCMCglmmRAM)
#making dataset ()


# need to compare the real 
mcmc_dat <- subset()
prior <- 
mcmc.thresh <- MCMCglmm(cbind() ~ trait -1, burnin =10000, nitt = 2e+05,
                     random = ~corg(trait):animal, rcov=~corg(trait):units, pedigree= td$phy, 
                     reduced = TRUE, data = td_glm,  prior = prior1, pr = TRUE, thin = 100,
                     pl = TRUE, family = c("categorical", "gaussian", "gaussian", "gaussian", "gaussian"))

y_dat <- as.factor(setNames(d_dat[,1], rownames(d_dat)))
dotTree(td$phy, y_dat, colors = c("blue","red"))

