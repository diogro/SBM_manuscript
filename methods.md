---
title: Gene co-expression and the Stochastic Block Model
bibliography: [./references.bib]
---

# Methods

## Gene clusterings

### Gene network

To model the structure of the gene co-expression correlation matrix, we use a network based approach using a Stochastic Block Model (SBM), which we describe below [@Peixoto2018-or]. This method attempts to uncover meaningful and statistically supported groups of nodes in a network based on their similarity with regards to their connections with other nodes. To do this, we start by generating a co-expression network graph, in which the genes are the nodes of the graph, and the inverse hyperbolic tangent transformed Spearman correlations between gene expressions are the weights in a fully connected graph. In theory, we could then procede using this full graph, but this is computationally too expensive. So, to reduce the size and connectivity of the co-expression graph, we impose a stringent false discovery rate (FDR) cut-off on the edges, removing edges with a large p-value associated with the correlation between the previously connected genes. As edges are removed, some genes with low or non-significant correlations become disconnected from the rest of the network and can be removed. By gradually reducing the FDR threshold, we reduce the size and density of the gene network, until we arrive at a viable set of genes and connections with which to fit the SBM. After some experimentation, we decided that an FDR of 1e-6 kept a large number of genes (~3000 for the body and ~4000 for the head) and allowed us to fit the SBM in days instead of weeks. This is a much larger number of meaningfully clustered genes than any other method (can’t remember the citations we had here…). Also, the weakly connected genes that we removed are not expected to significantly affect the partition of the network. (In other gene clustering methods, like WGCNA [@Langfelder2008-qa] or MMC [@Stone2009-hv], these weakly connected genes are usually placed in the "other" group and not used in further analysis)

### Stochastic Block Model

The Stochastic Block Model is a Bayesian generative model that attempts to find a partition of the gene network that minimizes the description length of the network. Broadly speaking, this is achieved by dividing the network into groups of genes, called blocks, and modeling the weight and existence of a link between two genes in a network solely on their belonging to a particular block. So, genes with similar patterns of connections tend to be clustered in the same block. The block structure that minimizes the description length of the network is equivalent to maximizing the posterior probability of the partition under the SBM. If $b$ is the modular partition, $A$ the weigthed gene network, we have: 

$$
P(b|A, \theta) = \frac{P(A|b, \theta) P(b)}{P(A)}
$$

### Nested SBM

The nested SBM uses a series of non-parametric hierarchical priors that greatly increase the resolution of block partition. This nested structure allows for the identification of more and smaller blocks that are statistically supported than other clustering methods [@Peixoto2017-zw]. This is achieved by treating the gene block partition as the nodes in a nested series of networks, which are then clustered using the same method. So, the genes are clustered in blocks, and these blocks are also clustered in a higher level blocks, and so on, as required to minimize the description lenght of the gene network. The model estimates the number of levels in the hierarchy and the number of blocks in each level. Since the model is generative, we can use posterior samples of the partitions to quantify error in any quantity estimated by the model, like the number of levels in the hierarchy, or the number of blocks at each level. For details in the implementation of the SBM, see [@Peixoto2017-zw] and [@Peixoto2018-or]. All SBM were fit using the graph-tool python library [@peixoto_graph-tool_2014]. 

### Edge weights

The weights on the edges, estimated from the gene expression correlations, can be modeled in the SBM using several different distributions. When using correlation, which are continuous numbers that vary between -1 and 1, it is  natural to use some transformation to map the correlations onto the real numbers. To do this, we use an the arctanh transformed correlations as the edge weights and model these weights using Gaussian distributions. In the SBM, the weights are modeled in much the same way as the links between networks, in that the the mean and the variance of the observed edge weigths between two blocks are a function only of the block structure, i.e., genes in the same block have a similar probability of being connected to other genes and the value of the weights in these edges come from the same distribution.

### Blocks vs. modules

Several co-expression methods attempt find group structure by maximizing modularity [@Newman2006-fv], with the objective of clustering genes into assortative modules: groups of genes that are more associated with each other than with other genes. Clustering genes in thighly correlated modules aligns with our intuition that groups of genes performing similar functions should be highly correlated, but it is not necessary for a network to be organized in this assortative fashion. The SBM is different from other clustering methods in that it does not attempt to find assortative modules. Instead, any information brough by the structural similarity between genes in the network can potentially be used to inform the clustering. To be sure, the SBM can capture an assortative modular pattern if it is present, but it is general enough to also capture other network organizations [@Zhang2020-up].  Assortative module maximization also has other know practical problems, paradoxically being prone to both overfitting (finding modular community structure where there is none) and under-fitting (failing to find modular structure in large networks, due to a problem know as the resolution limit [@Fortunato2007-ao]). The SBM avoids both of these problems, only revealing blocks when they are statistically supported and being able to detect a large number of groups when it is warranted [@Zhang2020-up]. In the context of the SBM, assortivity is not the main driver of the clustering, but it can be used to interpret the partitioning we obtain from the clustering. We can measure the assortivity of each group, and the global modularity of the resulting partitioning, and use this information in our analysis. 

### GO enrichment

 We ran GO in all blocks at all levels of the block hierarchy using R packages ClusterProfiler and XGR...

# Results

## SBM

In both head and body we identified a nested partition with 4 levels, with 2 blocks at level 4 (the coarsest), 4 in level 3, 12 (head) and 13 (body) at level 2, and 56 (head) and 53 (body) at level 1. 


### GO

  Tissue     Level 1          Level 2             Level 3    Level 4
--------    -------------   -----------------  ----------- -----------
head         73\% (41/56)     100\% (12)        100\% (4)    100\% (2)
body         77\% (41/53)     92\%  (12/13)     100\% (4)    100\% (2)
      
Table:  Number of blocks at each level that show significant GO enrichment at the 5\% FDR level.

### WGCNA comparison (this will probably be a separate supplementary document)

We compare the gene blocks from the SBM to the modules infered from WGCNA. WGCNA identifies 10 modules in the body and 20 in the head. The "other" module in WGCNA is composed of low degree genes that are not thighly linked to anything, and so are not included in any of the modules. SBM has no difficulty clustering these genes (due to the Degree Correction, not sure we want to get into this). Some WGCNA modules are remenicent of he recovered blocks, but several other things are included in the same group, while the SBM separates them more finelly. Enrichment in the WGCNA modules is less specific, and of the 20 reco

# References





