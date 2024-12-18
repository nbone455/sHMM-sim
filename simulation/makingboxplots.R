# Script to visualize model comparison results from four simulation scenarios
# Creates violin plots showing AIC differences between models 
# Data comes from previously saved simulation results (heatmap1-4)

#------------------------------------------------------------------------------
# Load required packages
#------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

#------------------------------------------------------------------------------
# Helper function to process heatmap results into plottable format
#------------------------------------------------------------------------------
  # Takes heatmap matrix and converts to long format data frame
  # heat: Matrix containing AIC differences from simulations
  # Returns: Data frame with columns for AIC, sample size (ntaxa), and number of states
make_box <- function(heat) {
  h1 <- heat
    
    # Process 4-state results
  h1_4 <- lapply(h1, '[', ,1)
  h1_4 <- unlist(h1_4)
  h1_4 <- data.frame(h1_4)
  h1_4$ntaxa <- "NA"
    # Add sample size information (50,100,200 taxa)
  h1_4$ntaxa[seq(from = 1, to = length(h1_4$ntaxa), by = 3)] <- 50
  h1_4$ntaxa[seq(from = 2, to = length(h1_4$ntaxa), by = 3)] <- 100
  h1_4$ntaxa[seq(from = 3, to = length(h1_4$ntaxa), by = 3)] <- 200
  colnames(h1_4) <- c("aic", "ntaxa")
  h1_4$nstate <- 4
  
  ggplot(h1_4, aes(x = ntaxa, y = aic)) +
    geom_boxplot()
  
  # 6 state 
  
  h1_6 <- lapply(h1, '[', ,2)
  h1_6 <- unlist(h1_6)
  
  h1_6 <- data.frame(h1_6)
  h1_6$ntaxa <- "NA"
  h1_6$ntaxa[seq(from = 1, to = length(h1_6$ntaxa), by = 3)] <- 50
  h1_6$ntaxa[seq(from = 2, to = length(h1_6$ntaxa), by = 3)] <- 100
  h1_6$ntaxa[seq(from = 3, to = length(h1_6$ntaxa), by = 3)] <- 200
  colnames(h1_6) <- c("aic", "ntaxa")
  h1_6$nstate <- 6
  
  ggplot(h1_6, aes(x = ntaxa, y = aic)) +
    geom_boxplot()
  
  # 8 state 
  
  h1_8 <- lapply(h1, '[', ,3)
  h1_8 <- unlist(h1_8)
  
  h1_8 <- data.frame(h1_8)
  h1_8$ntaxa <- "NA"
  h1_8$ntaxa[seq(from = 1, to = length(h1_8$ntaxa), by = 3)] <- 50
  h1_8$ntaxa[seq(from = 2, to = length(h1_8$ntaxa), by = 3)] <- 100
  h1_8$ntaxa[seq(from = 3, to = length(h1_8$ntaxa), by = 3)] <- 200
  colnames(h1_8) <- c("aic", "ntaxa")
  h1_8$nstate <- 8
  
  ggplot(h1_8, aes(x = ntaxa, y = aic)) +
    geom_boxplot()
  
  # 10 state 
  
  h1_10 <- lapply(h1, '[', ,4)
  h1_10 <- unlist(h1_10)
  
  h1_10 <- data.frame(h1_10)
  h1_10$ntaxa <- "NA"
  h1_10$ntaxa[seq(from = 1, to = length(h1_10$ntaxa), by = 3)] <- 50
  h1_10$ntaxa[seq(from = 2, to = length(h1_10$ntaxa), by = 3)] <- 100
  h1_10$ntaxa[seq(from = 3, to = length(h1_10$ntaxa), by = 3)] <- 200
  colnames(h1_10) <- c("aic", "ntaxa")
  h1_10$nstate <- 10
  
  ggplot(h1_10, aes(x = ntaxa, y = aic)) +
    geom_boxplot()
  
  # 12 state 
  
  h1_12 <- lapply(h1, '[', ,4)
  h1_12 <- unlist(h1_12)
  
  h1_12 <- data.frame(h1_12)
  h1_12$ntaxa <- "NA"
  h1_12$ntaxa[seq(from = 1, to = length(h1_12$ntaxa), by = 3)] <- 50
  h1_12$ntaxa[seq(from = 2, to = length(h1_12$ntaxa), by = 3)] <- 100
  h1_12$ntaxa[seq(from = 3, to = length(h1_12$ntaxa), by = 3)] <- 200
  colnames(h1_12) <- c("aic", "ntaxa")
  h1_12$nstate <- 12
  
  ggplot(h1_12, aes(x = ntaxa, y = aic)) +
    geom_boxplot()
  
  # 14 state 
  
  h1_14 <- lapply(h1, '[', ,4)
  h1_14 <- unlist(h1_14)
  
  h1_14 <- data.frame(h1_14)
  h1_14$ntaxa <- "NA"
  h1_14$ntaxa[seq(from = 1, to = length(h1_14$ntaxa), by = 3)] <- 50
  h1_14$ntaxa[seq(from = 2, to = length(h1_14$ntaxa), by = 3)] <- 100
  h1_14$ntaxa[seq(from = 3, to = length(h1_14$ntaxa), by = 3)] <- 200
  colnames(h1_14) <- c("aic", "ntaxa")
  h1_14$nstate <- 14
  
  ggplot(h1_14, aes(x = ntaxa, y = aic, fill = nstate)) +
    geom_boxplot()
  
  
  # 16 state 
  
  h1_16 <- lapply(h1, '[', ,4)
  h1_16 <- unlist(h1_16)
  
  h1_16 <- data.frame(h1_16)
  h1_16$ntaxa <- "NA"
  h1_16$ntaxa[seq(from = 1, to = length(h1_16$ntaxa), by = 3)] <- 50
  h1_16$ntaxa[seq(from = 2, to = length(h1_16$ntaxa), by = 3)] <- 100
  h1_16$ntaxa[seq(from = 3, to = length(h1_16$ntaxa), by = 3)] <- 200
  colnames(h1_16) <- c("aic", "ntaxa")
  h1_16$nstate <- 16
  
  ggplot(h1_16, aes(x = ntaxa, y = aic, fill = nstate)) +
    geom_boxplot()
  
  
  h1_full <- rbind(h1_4, h1_6, h1_8,h1_10, h1_12,h1_14, h1_16)
  
  
  # heatmap 1, 0 aic would imply no diff between threshold and
  # mk, nbegative is preferred threshold model and posi is preferred
  # Mk model
  return(h1_full)
}

library(ggplot2)
library(dplyr)

setwd("~/Documents/WORK/Paper1/")
h1 <- readRDS("./boxplots_sims/heatmap_results/heatmap1_aic_results.rds")

h1_full <- make_box(h1)

h1_full$ntaxa <- factor(h1_full$ntaxa , levels=c("50", "100", "200"))
h1_fullb <- mutate(h1_full, sig=abs(h1_full$aic)>4)

# all of them
ggplot(h1_fullb, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  facet_wrap(~nstate, scale="fixed") + theme(legend.position="none")


# individual 
## just four 
h1_four <- filter(h1_fullb, nstate == 4)
four <- ggplot(h1_four, aes(x = ntaxa, y = aic)) +
  ylim(-20,20) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = four, filename = "h1_four.png",dpi=900)

## just six
h1_six <- filter(h1_fullb, nstate == 6)
six <- ggplot(h1_six, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = six, filename = "h1_six.png",dpi=900)

## just eight 
h1_eight <- filter(h1_fullb, nstate == 8)
eight <- ggplot(h1_eight, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = eight, filename = "h1_eight.png",dpi=900)

## just ten
h1_ten <- filter(h1_fullb, nstate == 10)
ten <- ggplot(h1_ten, aes(x = ntaxa, y = aic)) +
  ylim(-20,20) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = ten, filename = "h1_ten.png",dpi=900)

## just 12
h1_twelve <- filter(h1_fullb, nstate == 12)
twelve <- ggplot(h1_twelve, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = twelve, filename = "h1_twelve.png",dpi=900)

## just 14 states
h1_fourteen <- filter(h1_fullb, nstate == 14)
fourteen <- ggplot(h1_fourteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  ylim(-20,20)+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = fourteen, filename = "h1_fourteen.png",dpi=900)

## just 16

h1_sixteen <- filter(h1_fullb, nstate == 16)
sixteen <- ggplot(h1_sixteen, aes(x = ntaxa, y = aic)) +
  ylim(-20,20) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = sixteen, filename = "h1_sixteen.png",dpi=900)

# H2

h2 <- readRDS("./boxplots_sims/heatmap_results/heatmap2_aic_results.rds")
h2_full <- make_box(h2)  
h2_full$ntaxa <- factor(h2_full$ntaxa , levels=c("50", "100", "200"))
h2_fullb <- mutate(h2_full, sig=abs(h2_full$aic)>4)

ggplot(h2_fullb, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  facet_wrap(~nstate, scale="fixed") + theme(legend.position="none")

# individual 
## just four 
h2_four <- filter(h2_fullb, nstate == 4)
four <- ggplot(h2_four, aes(x = ntaxa, y = aic)) +
  ylim(-20,20) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = four, filename = "h2_four.png",dpi=900)

## just six
h2_six <- filter(h2_fullb, nstate == 6)
six <- ggplot(h2_six, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = six, filename = "h2_six.png",dpi=900)

## just eight 

h2_eight <- filter(h2_fullb, nstate == 8)
eight <- ggplot(h2_eight, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  ylim(-20,20) +
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = eight, filename = "h2_eight.png",dpi=900)

## just ten
h2_ten <- filter(h2_fullb, nstate == 10)
ten <- ggplot(h2_ten, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = ten, filename = "h2_ten.png",dpi=900)

## just 12

h2_twelve <- filter(h2_fullb, nstate == 12)
twelve <- ggplot(h2_twelve, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = twelve, filename = "h2_twelve.png",dpi=900)

## just 14

h2_fourteen <- filter(h2_fullb, nstate == 14)
fourteen <- ggplot(h2_fourteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = fourteen, filename = "h2_fourteen.png",dpi=900)

## just 16

h2_sixteen <- filter(h2_fullb, nstate == 16)
sixteen <- ggplot(h2_sixteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = sixteen, filename = "h2_sixteen.png",dpi=900)


# H3


h3 <- readRDS("./boxplots_sims/heatmap_results/heatmap3_aic_results.rds")
h3_full <- make_box(h3)
h3_fullb <- mutate(h3_full, sig=abs(h3_full$aic)>4)
h3_alt <- h3_fullb
h3_alt$over <- 0

# fill in the 'over' column with a 1 if the absolute value 
# of the 'aic' column is more than 10000, or a 0 if less than
h3_alt$over<- case_when(
  abs(h3_alt$aic) > 10000 ~ 1,
  abs(h3_alt$aic) < 10000 ~ 0
)

h3_alt <- h3_alt[h3_alt$over == 0, ] 

# individual 
## just four 
h3_four <- filter(h3_fullb, nstate == 4)
four <- ggplot(h3_four, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = four, filename = "h3_four.png",dpi=900)

## just six
h3_six <- filter(h3_fullb, nstate == 6)
six <- ggplot(h3_six, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = six, filename = "h3_six.png",dpi=900)

## just eight 

h3_eight <- filter(h3_fullb, nstate == 8)
eight <- ggplot(h3_eight, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  ylim(-20,20) +
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = eight, filename = "h3_eight.png",dpi=900)

## just ten
h3_ten <- filter(h3_fullb, nstate == 10)
ten <- ggplot(h3_ten, aes(x = ntaxa, y = aic)) +
  ylim(-30,30)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = ten, filename = "h3_ten.png",dpi=900)

## just 12

h3_twelve <- filter(h3_fullb, nstate == 12)
twelve <- ggplot(h3_twelve, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = twelve, filename = "h3_twelve.png",dpi=900)

## just 14

h3_fourteen <- filter(h3_fullb, nstate == 14)
fourteen <- ggplot(h3_fourteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = fourteen, filename = "h3_fourteen.png",dpi=900)

## just 16

h3_sixteen <- filter(h3_fullb, nstate == 16)
sixteen <- ggplot(h3_sixteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = sixteen, filename = "h3_sixteen.png",dpi=900)



## need to drop crazy outliers
ggplot(h3_alt, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig), binaxis='y', stackdir='center', dotsize=0.8) +
  facet_wrap(~nstate, scale="fixed") + theme(legend.position="none")






h4 <- readRDS("./boxplots_sims/heatmap_results/heatmap4_aic_results.rds")
h4_full <- make_box(h4)
h4_fullb <- mutate(h4_full, sig=abs(h4_full$aic)>4)


# individual 
## just four 
h4_four <- filter(h4_fullb, nstate == 4)
four <- ggplot(h4_four, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = four, filename = "h4_four.png",dpi=900)

## just six
h4_six <- filter(h4_fullb, nstate == 6)
six <- ggplot(h4_six, aes(x = ntaxa, y = aic)) +
  ylim(-20,20)+
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = six, filename = "h4_six.png",dpi=900)

## just eight 

h4_eight <- filter(h4_fullb, nstate == 8)
eight <- ggplot(h4_eight, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = eight, filename = "h4_eight.png",dpi=900)

## just ten
h4_ten <- filter(h4_fullb, nstate == 10)
ten <- ggplot(h4_ten, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = ten, filename = "h4_ten.png",dpi=900)

## just 12

h4_twelve <- filter(h4_fullb, nstate == 12)
twelve <- ggplot(h4_twelve, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = twelve, filename = "h4_twelve.png",dpi=900)

## just 14

h4_fourteen <- filter(h4_fullb, nstate == 14)
fourteen <- ggplot(h4_fourteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = fourteen, filename = "h4_fourteen.png",dpi=900)

## just 16

h4_sixteen <- filter(h4_fullb, nstate == 16)
sixteen <- ggplot(h4_sixteen, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(device = 'png', plot = sixteen, filename = "h4_sixteen.png",dpi=900)


ggplot(h4_fullb, aes(x = ntaxa, y = aic)) +
  geom_violin(aes())+
  geom_dotplot(aes(fill=sig, color=sig, ), binaxis='y', stackdir='center', dotsize=0.8) +
  facet_wrap(~nstate, scale="fixed") + theme(legend.position="none")

