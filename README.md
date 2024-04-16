# sHMM-sim
Repository for reproducing all the analyses of the "detecting threshold traits from discrete data using Structured hidden Markov models" paper 

**Abstract**: Testing hypotheses of trait evolution often requires the simplification of complex phenotypes with continuous variation into discrete categories or states. While this approach is sufficient for many traits and evolutionary questions, it can also conceal the underlying evolutionary patterns and introduce subjective biases into statistical analyses. This study explores the efficacy of structured hidden Markov models (SHMMs) in testing hypotheses related to the underlying complexity of those simplified categories. Traditional methods can inadvertently introduce assumptions about a trait's underlying variation that conflict with its phylogenetic history, potentially diminishing the power of inferences. SHMMs, by incorporating unobserved 'hidden' states, offer a nuanced approach to model the evolution of discrete characters. Through a series of simulations, we demonstrate that SHMMs, when equipped with an adequate number of hidden states and a sufficiently large phylogeny, outperform standard Markov models in identifying true threshold traits--discretized versions of traits with one-dimensional continuous variation. An empirical examination of primate molar evolution further underscores the potential of SHMMs in deciphering the intricate dynamics of trait evolution. Our findings advocate for a more structured approach to character state modeling and hypothesis testing, emphasizing the importance of validating expert knowledge in trait categorization

# Simulations

four scenarios are set up: 

## heatmap 1: true model threshold 
normal Markov equal rates vs threshold-approximation hidden state model equal rates

## heatmap 2: true model Mk 
normal Markov equal rates vs threshold-approximation hidden state model equal rates

## heatmap 3: true model threshold
generalized hidden Markov model vs threshold-approximation hidden state model equal rates

## heatmap 4: true model is a "almost" a threshold model
generalized hidden Markov model vs threshold-approximation hidden state model equal rates
