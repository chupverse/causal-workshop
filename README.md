# Best practices and recent advances in causal analyses

This page is devoted to the practical phase of the Inserm workshop related to "Best practices and recent advances in causal analyses". Four session are organised. This repsitory aims to share the supports used during the three training days.

### Variable selection with DAGitty (repository [session1](https://github.com/chupverse/causal-workshop/tree/main/session1))

In observational epidemiological studies, it is important to identify confounding factors to estimate the causal effect of exposure on the health event. There are various methods for selecting these factors, including selection based on graphs representing causal hypotheses between the different factors, exposure, and the health event (hypotheses developed from the literature and knowledge of the context). This identification strategy is implemented in the free software DAGitty. In this session, we will use concrete examples from epidemiological studies to illustrate the construction of a causal graph and its implementation in DAGitty.

### G-computation and propensity score weighting  (repository [session2](https://github.com/chupverse/causal-workshop/tree/main/session2))

Weighting consists of individually predicting the probability of receiving treatment and then decreasing/increasing the contributions of overrepresented/underrepresented individuals. The G-computation consists of individually predicting the responses to treatment and then estimating the marginal effect by averaging. We will first explain in detail the two methods. Then, we will use the methods on the same application to compare the pros and cons.

### Mediation analysis with the R CMAverse package  (repository [session3](https://github.com/chupverse/causal-workshop/tree/main/session3))

Mediation analysis breaks down the total effect of exposure on a variable of interest into its direct effect and its effect mediated by one or more intermediate variables. We will discover the CMAverse R package, which allows mediation analyses to be performed in various contexts (different types of effects, causal diagrams and types of variables). Using concrete examples, we will discuss the assumptions underlying the identification of causal effects, implement their estimation, interpret the causal contrasts obtained, and discuss the results.

### Mediation analyses with the ltmle and medoutcon packages   (repository [session4](https://github.com/chupverse/causal-workshop/tree/main/session4))

The ltmle package allows you to estimate the causal effect of repeated exposures over time with three possible estimators: g-computation by iterative expectations, weighting, and Targeted Maximum Likelihood Estimation (TMLE). In the context of mediation analyses, this package can be used to estimate controlled direct effects, even in the presence of time-dependent confounding influenced by initial exposure (recanting witness). We will review the principles of double robustness estimation using TMLE and work through a simple illustrative example to understand how this package is implemented. When our objective is to estimate a natural direct or indirect effect rather than a controlled direct effect, there are fewer tools available for obtaining a double robust estimate. We will use the medoutcon package, which provides a one-step estimator and a TMLE-based estimator.
