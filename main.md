---
title: Clustering using the Stochastic Block Model shows that gene expression networks are not assortative
author:
- Diogo Melo
- Luisa Pallares
- Julien Ayroles
output: pdf_document
geometry:
- top=20mm
- left=25mm
- right=25mm
- bottom=20mm
header-includes:
- \usepackage[backref=true,style=authoryear]{biblatex}
- \DefineBibliographyStrings{english}{backrefpage = {page}, backrefpages = {pages}}
- \usepackage{multicol}
- \usepackage{setspace}
- \usepackage{float}
- \usepackage{afterpage}
- \usepackage{stfloats}
- \usepackage{graphicx}
- \newcommand{\hideFromPandoc}[1]{#1}
- \hideFromPandoc{ \let\Begin\begin \let\End\end}
link-citations: yes
mainfont: Skolar PE TEST Regular
mainfontoptions:
- Numbers=Lowercase
- Numbers=Proportional
csl: ./evolution.csl
sansfont: Skolar Sans PE TEST
bibliography: ./references.bib
---

# Intro

Gene co-expression networks allow for the associations between genes to inform our biological understanding of cell and organismal function. Associations between gene expressions can be indicative of common function, and the number of connections can be an indicator of central or regulatory genes [@Van_Dam2018-nf]. Due to the large dimensionality of gene expression data, often composed of several thousands of gene expressions, one major tool in the analysis of co-expression is gene clustering: separating the genes into related groups, which can then be explored separately [@Dhaeseleer2005-jv]. This drastically reduces the number of genes we need to consider at the same time, and allows for the identification of hubs or centrally connected genes that can then be used to inform further experimental validation [@Langfelder2008-qa; @Imenez_Silva2017-ic]. 

The question is, given a co-expression network, how should we cluster the genes? The general idea behind several methods is to look for similar genes, as these are expected to be involved in related biological functions, but several definitions of similarity have been used. The most basic measure of similarity borrows from classical morphological integration theory, and attempts find gene modules based on their correlations. In this context, genes in the same modules are expected to be highly correlated and perform similar functions, while genes in different modules are expected to have low correlations [@Olson1958-fh; @Magwene2001-an; @Wagner2007-jt]. Here, we refer to this classic pattern of groups of genes that are more associated with each other than with other genes as assortativity, and to the groups as assortative modules. Other methods use the correlations to create other measures of similarity, which are then used as inputs to clustering algorithms. Weighted correlation network analysis [WGCNA, @Langfelder2008-qa] uses a power transformation of the correlations between genes expressions (or a topological similarity measure built with these transformed correlations) [@Zhang2005-kh; @Dong2007-ff] as a similarity measure that is then separated into assortative modules using hierarchical clustering. One of the main objectives of WGCNA is finding hub genes, which have high connectivity withing modules and are clearly identified by the hierarchical clustering. Other methods borrow from network analysis and attempt to explicitly maximize the Newman Modularity [@Newman2006-fv] of the weighted gene network. Modulated Modularity Clustering [MMC, @Stone2009-hv] uses an adaptive algorithm to find a non-linear distance between genes based on their correlations that maximizes the number of modules. All these methods impose an assortative structure on the gene expression network, in which similar genes are expected to be correlated with each other. 

These clustering approaches come with important downsides. Using WGCNA involves manually tuning several parameters, and for this tuning the method leans heavily on the expectation that gene co-expression networks should be approximately scale-free [@Dong2007-ff; @Bergmann2004-vw; @Jeong2000-xe], but, despite its popularity, this expectation might be somewhat unwarranted [@Khanin2006-ve; @Stumpf2005-ed]. Even with optimal parameters, WGCNA fails to assign a substantial proportion of genes to any module. While WGCNA is efficient in finding hub genes, if some functional gene group does not have a hub or has low average similarities, this group will never be identified. Both limitations potentially leave biological insight on the table by ignoring network structures that are different from what the method expects. Methods that use modularity maximization, like MMC, are subject to know statistical problems, surprisingly being prone to both overfitting (finding modular community structure where there is none [@Guimera2004-jq]) and under-fitting (failing to find modular structure in large networks, due to a problem know as the resolution limit [@Fortunato2007-ao]). 

The co-expression methods described above attempt find group structure by clustering genes into assortative modules. Clustering genes in tightly correlated modules aligns with our intuition that groups of genes performing similar functions should be highly correlated, but it is not necessary for a network to be organized in this assortative fashion, and so using methods that assume assortativity will necessarily ignore alternative organizations. Here, we use a more general measure of similarity can be used to find meaningful genes groups that are not necessarily assortative, with high correlation between groups and biological function. To achieve this, we use a weighted nested degree corrected stochastic block model (wnDC-SBM, or SBM for brevity) [@Peixoto2017-zw; @Peixoto2018-or], which has shown promising results in similar applications [@Baum2019-ty; @Morelli2020-ge]. The SBM is different from other clustering methods in that it does not attempt to find assortative modules. Instead, any information contained in the structural similarity between genes in the network can potentially be used to inform the clustering. To be sure, the SBM can capture an assortative modular pattern if it is present, but it is general enough to also capture other network organizations [@Zhang2020-up]. Even if, in the context of the SBM, assortativity is not the main driver of the gene partitioning, it can still be used to interpret the groups we obtain. We show that the SBM, a model with no free parameters, can find many more meaningful groups than competing methods, as revealed by highly specific gene ontology enrichment, and show that the expected assortativity of biologically meaningful groups is not necessarily present.


# Methods

## Gene expression measures

We measure gene expression using RNASeq in a large, outbred population of D. Melanogaster. We measure expression separately for the body and the head. 

## Gene selection

Using the gene expression measures for both tissues we generate co-expression network graphs. In theory, we could proceed using a full network in which all pairs of genes are connected but fitting the SBM with this fully connected graph is computationally too expensive. So, to reduce the size and connectivity of the network, we impose a stringent Benjamini-Hochberg false discovery rate (FDR) cut-off on the edges, removing edges with a large p-value associated with the correlation between the corresponding genes. As edges are removed, some genes with only non-significant correlations become disconnected from the rest of the network and can be removed. By gradually reducing the FDR threshold, we reduce the size and density of the gene network, until we arrive at a viable set of genes and connections with which to fit the SBM. After some experimentation, we decided that an FDR of 1e-2 for the head and 1e-3 for the body kept most of the genes (head:5261, body:5124) and reduced the graph density to a manageable level. This set of genes is used in all the different clustering methods. 

## Edge weights

Each method uses different edge weights for the network graph. Both WGCNA and MMC can use the fully connected graph, so we maintain all edges in these methods. We use the TOM similarity in WGCNA, and the Spearman correlation derived distance in MMC. For the SBM, we use the low-density graph described above, with the edge weights given by the inverse hyperbolic tangent transformed Spearman correlations between gene expressions. This transformation allows the edge weights to be modeled by normal distributions in the SBM, as we discuss below.

## Stochastic Block Model

The Weighted Nested Degree Corrected Stochastic Block Model [@Karrer2011-vp; @Peixoto2017-zw] is a Bayesian generative model that attempts to find the partition with the highest posterior probability given the observed network and edge weights. Broadly speaking, this is achieved by dividing the network into groups of genes, called blocks, and modeling the weight and existence of a link between two genes in a network solely on their belonging to a particular block. So, genes with similar patterns of connections tend to be clustered in the same block. The degree correction refers to a modification of the standard Stochastic Block Model that allow genes with different connectivity to be clustered in the same block [see @Peixoto2017-zw for details].

If $b$ is a particular partition of the genes in the weighted gene network $A$, we write a model that generates $A$ with probability given by $P(A| b, \theta)$, where $\theta$ stands in for any extra parameter we need besides the group partition $b$. With this model, we can write the posterior probability of the block partition $b$ given the observed network:

$$
P(b|A) = \frac{P(A|\theta, b)P(\theta, b)}{P(A)}
$$

where P(A) is a normalization constant. As for the additional parameters $\theta$, the formulation used here, from @Peixoto2018-or, uses hard constraints such that there is only one choice of $\theta$ that is compatible with $A$ and $b$, which means that the model has no free parameters. We can then search for the partition $b$ that maximizes $P(b|A)$ using computational methods, like Markov Chain Monte Carlo (MCMC) methods.

### Description length

The posterior probability of the block partition can be written as:

$$
P(b|A) \propto \exp(-\Sigma)
$$

Where $\Sigma = -log[P(A|\theta, b)] - log[P(\theta,b)]$ is called the description length of the gene network $A$, and has an information theoretic interpretation, being the amount of information required to encode the network given $\theta$ and $b$. So finding the partition that maximizes the posterior probability is the same as minimizing the description length, or, in other words, the chosen partition $b$ allows us to describe the network using less information. 

The two terms in $\Sigma$ also allow us to understand why this method offers an intrinsic protection against overfitting. The first term $log[P(A|\theta, b)]$ corresponds to the log likelihood of the observed network. Increasing the number of blocks allows this likelihood to increase as the degrees of freedom of the model increase. But, the second term, $log[P(\theta,b)]$ functions as a penalty that decreases for complex models with many blocks, and the description length cannot decrease for overly complex models that have more blocks than warranted by the data. So, the selected partition with the minimum description length will necessarily be the simplest partition with similar explanatory power, avoiding overfitting and fully based on the available statistical evidence. For example, he SBM would not detect modules that appear in random networks due to statistical fluctuations [@Guimera2004-jq; @Zhang2020-up]. We can also use the description length as a principled method for comparing models that simultaneously considers fit to data and model complexity.

### Weighted SBM

The weights on the edges can be modeled in the SBM using several different distributions. When using correlation, which are continuous numbers that vary between -1 and 1, it is natural to use some transformation to map the correlations onto the real numbers. To do this, we use arctanh transformed correlations as the edge weights and model these weights using normal distributions. In the SBM, the weights are modeled in much the same way as the links between networks, in that the the mean and the variance of the observed edge weights between two blocks are a function only of the block structure, i.e., genes in the same block have a similar probability of being connected to other genes and the value of the weights in these edges come from the same distribution.

### Nested SBM

The nested SBM uses a series of non-parametric hierarchical priors that greatly increase the resolution of block partition. This nested structure allows for the identification of more and smaller blocks that are statistically supported than other clustering methods [@Peixoto2017-zw]. This is achieved by treating the gene block partition as the nodes in a nested series of networks, which are then clustered using the same method. So, the genes are clustered in blocks, and these blocks are also clustered in a higher-level blocks, and so on, as required to minimize the description length of the gene network. The model estimates the number of levels in the hierarchy and the number of blocks in each level. Since the model is generative, we can use posterior samples of the partitions to quantify error in any quantity estimated by the model, like the number of levels in the hierarchy, or the number of blocks at each level. For details in the implementation of the SBM, see [@Peixoto2017-zw] and [@Peixoto2018-or]. All SBM were fit using the graph-tool python library [@peixoto_graph-tool_2014]. 

### Modularity and Assortativity

Instead of attempting to maximize modularity, when using the nested SBM we can ask if the inferred partition is modular or not by calculating the Newman modularity at each level of the nested hierarchy. Modularity is calculated at each nested level using:

$$
M = \frac{1}{2E} \sum_r e_{rr} - \frac{e_{r}^2}{2E}
$$

where $e_{rs}$ is the sum of edge weights between groups $r$ and $s$, $e_{r} = \sum_s  e_{rs}$ and $E$ is the sum of all weights. We further decompose the contribution of each Level-1 block to the modularity, by defining the assortativity of a block as:

$$
q_r = \frac{B}{2E} \left ( e_{rr} - \frac{e_{r}^2}{2E} \right )
$$

where $B$ is the number of blocks. Using this definition, $M = \frac{1}{B} \sum_r  q_r$, and modules with positive assortativity contribute to increasing modularity, while blocks with negative assortativity decrease it. Assortativity values vary between -1 for fully disassortative block and 1 for a fully assortative one.

## WGCNA and MMC

We use WGCNA to cluster the genes into modules using the topological overlap measure (TOM) similarity with a soft threshold of 6 in a signed similarity measure. WGCNA produces modules by cutting the hierarchical clustering tree at different heights, and we use the dynamic cutting option to create the modules. We use a signed network (as opposed to ignoring the sign of the correlation between genes) because inspection of the gene network graph reveals large groups of genes linked by negative correlations in our data, suggesting large scale structure that would be obscured by using the unsigned method. Unsigned similarity has been shown to lead to more robust modules [@Mason2009-ej], and in tuning WGCNA we were able to cluster more genes and find more modules using the signed method. MMC has no option to use the sign of the correlation, so we use the absolute value of the Spearman correlations.

## Gene Ontology enrichment

We assess the biological relevance of the clustering obtained by each method by comparing their gene ontology (GO) enrichment. We filter enrichment using a Benjamini-Hochberg FDR rate of 5\%, with a minimum of 4 genes in the enriched set. We ran GO in all identified clusters using the clusterProfiler R package v4.2.2 [@Yu2012-tz].

# Results

## Clustering

Gene clustering for all methods is presented in table S1. Using the SBM, in both head and body we identified a nested partition with 5 levels, with 2 blocks at level 5 (the coarsest), 3 blocks at level 4, 6 (head) and 9 (body) in level 3, 21 (head and body) at level 2, and 82 (head) and 78 (body) at level 1. The block structure inferred by the SBM is shown in fig. @fig:Emats. In what follows, when discussing specific SBM blocks, we either explicitly define which level of the nested hierarchy we are referring to, or give the full path to a given block. So, level-1 block 12 in the head can also be referred as 12-7-2-2-1, and level-2 block 10 in the body is also 10-1-2-1. 

![SBM Level 1 blocks colored by their weighted connectivity, within and between blocks. Upper levels of the nested hierarchy are shown by the red lines.](figures/SBM_Ematrix.png){#fig:Emats}

WGCNA partitioned 2118 genes into 7 modules in the body and 1600 genes into 7 in the head. WGCNA did not cluster 3006 genes in the body and 3661 in the head. Given that the number of modules in WGCNA is somewhat similar to the number of blocks at level 3 of the SBM, we compare these two partitions in fig. @fig:wgcna_compare. Overall, the partitions are different, but there are some common patterns. For example, Level-3 blocks 0, 2, 5, and 6 in the body are split between modules 3 and 4, and these blocks are all in the same Level-4 block 0, suggesting some similarity that could explain the WGCNA clustering. Blocks 7 and 9 are both fully assigned to module 2. Also in the body genes, we find a similar patter for Level-3 blocks 1, 3, and 4, which are mostly split between modules 1 and 2. In the head, Level-3 block 4 is all assigned to modules 1 and 2. Level-3 blocks 1 and 2 are split between modules 1 and 2, and both are in Level-4 block 2. So, while the clustering is different, WGCNA and the SBM do capture some common signal.

MMC failed to cluster most genes, placing almost all genes into the same large module. Given this poor performance on our data, we do not discuss MMC further and instead focus on comparing SBM to WGCNA.

![Comparison of the clustering in WGCNA and in levels 2 and 3 of the SBM for the gene expressions in the body (left) and the head (right). Each point corresponds to a gene. The x axis corresponds to the Level 2 SBM blocks, y axis the WGCNA modules. Colors correspond to the (coarser) level 3 of the SBM.](figures/WGCNA_comparison.png){#fig:wgcna_compare}

## Modularity and assortativity

Modularity and assortativity are markedly lower in the body (fig. @fig:modularity). Several blocks in the body have negative assortativity, and the maximum value of modularity is 0.035 at level 4 of the nested hierarchy. Even so, several blocks show GO enrichment across the distribution of assortativity. In the head, modularity is overall higher, with a peak at 0.14 in level 3. This is still a relatively low value, and illustrates how assuming the gene network should be modular can prevent us from finding an informative clustering. All but 5 blocks in the head show positive assortativity, and again GO enrichment is present across the range. 


![Assortativity in the SBM level-1 blocks and modularity for all nested levels. GO enriched blocks are show in yellow and appear all along the distribution of assortativity. Modularity is much higher in the head, and it peaks at level 3, dropping in upper levels. Body has a much higher number of non-assortative blocks and lower modularity at all levels. Modularity peaks at level 4 in the body, and drops strongly at level 5. Interestingly, the 4 most assortative blocks in the body are do not show significant GO enrichment.](figures/assortativity.png){#fig:modularity}

## Gene Ontology enrichment

The majority of blocks in SBM show some level of GO enrichment (Table 1). In particular, several of the Level-2 blocks show a remarkable consistency in their enrichment. For example, Level-3 block 0 in the head is related to neural signaling and sensory perception, with its daughter blocks at Level-2 (4 and 6) showing enrichment for: (Level-2 4) G protein-coupled receptor signaling pathway, detection of light stimulus, phototransduction (Level-2 6) synapse organization, axon development,  cell-cell signaling, behavior. Several of these enrichments are exclusive to one of the level-1 blocks. Perhaps the most surprising enrichment is for the blocks associated with translation, in which all of the ribosomal proteins are clustered almost exclusively in level-1 blocks in the body and in the head. We discuss some notable blocks below. Several other Level-2 and Level-1 blocks are readily identifiable as related to development, DNA transcription, cell respiration, cell cycle regulation, immune response, sugar metabolism, among others. All WGCNA modules show GO enrichment (but modules 5, 6 and 7 in the body show only 2 or less enriched terms). Supporting Information table 1 shows GO enrichment for all SBM blocks and WGCNA modules. 

### Notable individual clusters

Level 2 block 0-0-0 in the head is one of the easiest to interpret, being entirely related to brain tissue function. Figure @fig:go_map shows the top 10 GO categories for each of the level-1 blocks in block 0-0-0, and the most neuronal enriched WGCNA grouping, module 4. The SBM blocks separate vesicle exocytosis, neuronal differentiation, phototransduction, synaptic signaling, and, interestingly, there is a block related to alternative mRNA splicing, which is though to be more common in brain tissues [@Su2018-nz]. WGCNA module 4 recovers some of this enrichment, but in a less granular way. Only the phototransduction part of the enrichment is separated in module 6. Several of these level-1 blocks (especially 0, 7 and 24) are among the most assortative (see fig. @fig:modularity and following section), and so are prime candidates for detection in WGCNA. The alternative splicing module has a much lower assortativity, and so its not surprising that WGCNA could not detect it.

Some of the most specific enrichments in the SBM are the translation related blocks. In both body and head, ribosomal proteins are clustered in small and highly enriched level-1 blocks. Six level-1 blocks in the head and seven in the body are composed of virtually only ribosomal related protein genes (We show a subset of these in Figure @fig:go_translation). All are very specific, being composed of between 10-22 genes, and have low assortativity (mean: $0.007$, max: $0.07$, min: $-0.03$). There is no equivalent module in WGCNA, but all of the translation related genes are in the same much larger modules (module 1 in the head, 449 genes; and module 5 in the body, 324 genes), both of which show enrichment for translation but also for several other categories. The level-2 clustering of level-1 blocks in the SBM is also informative. In the head, all of the translation level-1 blocks are in their own level-2 blocks (1-1-1 and 6-1-1). In contrast, in the body, the level-1 translation block share a level-2 blocks with cell respiration blocks. In the body, blocks 1-1-1 and 11-1-1 are about evenly split into translation and mitochondrial respiration level-1 blocks (fig. @fig:body_11 shows one of these mixed level-2 blocks). WGCNA also places cell respiration related genes in the body on the same module 5.

The SBM is also able to capture exceedingly faint biological signals. For example, level-1 blocks 55 in the head and 31 in the body are similar in that they are visibly less connected to the rest of the blocks, as we can see in fig. @fig:Emats. They are both very small, being composed of 17 and 10 genes, respectively. This similarity in confirmed by the common GO enrichment for both blocks (immune response) and some common genes (Dpt - diptericin, Dro - drosocin).

  Tissue     Level 1          Level 2           Level 3        Level 4     Level 5
--------    --------------   ---------------  ------------- ------------- -------------
head         63\% (52/82)     95\% (20/21)     100\% (6/6)   100\% (3/3)   100\% (2/2)
body         68\% (56/78)     95\% (20/21)     100\% (9/9)   100\% (3/3)   100\% (2/2)
      
Table: Number of blocks at each level that show significant GO enrichment at the 5\% FDR level with a minimum of 4 genes in the enriched set.

![Enriched GO categories in a level-2 block in the head (0-0-0), related to neural signaling. Panels show the corresponding level-1 blocks. Labeled nodes are the top 10 GO categories by number of genes, gray unlabeled nodes are genes linked to the GO categories. The last panel shows the most similar WGCNA module, which also contains signaling related genes, but at a lower resolution and fails to cluster the phototransduction genes, which are in module 6.](figures/000_go_map.png){#fig:go_map}

![Enriched GO categories for a subset of the translation related level-1 blocks for the Head and Body. Labeled nodes are the top 10 GO categories by number of genes, gray unlabeled nodes are genes linked to the GO categories.](figures/Translation_go_map.png){#fig:go_translation}

![Enriched GO categories in a level-2 block in the body (11-1-1), related to translation and cell respiration. Panels show the corresponding level-1 blocks. Labeled nodes are the top 10 GO categories by number of genes, gray unlabeled nodes are genes linked to the GO categories.](figures/plot_11_1_1.png){#fig:body_11}


\newpage

# Discussion

Traits in an organism need to have some level of integration, of interdependence, in order to form a functioning individual. This necessary interaction between parts poses an important problem for understanding the evolution of complex traits, as inter-dependencies are expected to lead to important evolutionary restrictions [@Orr2000-gn]. Modularity, relative independence  between groups of complex traits, provides a simple solution to this problem as it allows organisms to maintain their function unchanged by coordinating simultaneous evolutionary changes in all related traits, while keeping unrelated traits undisturbed [@Ancel2000-vt; @Cheverud1996-jw; @Wagner2011-zx; @Wagner1996-ui]. The explanatory power of modularity as a unifying concept in several levels of organization has lead to its use in several different fields with great success [@Melo2016-yw; @Zelditch2021-ue], and there is no doubt that there are several biological systems that are indeed modular. This general usefulness has also informed much of our thinking on how complex traits should be structured, producing a large literature dedicated to finding modules and testing for their existence. However, modularity is not a necessary feature of biological organization [even in the case of evolvability, see @Pavlicev2011-xm], and restricting our attention to modular systems can blind us to alternative organizations. Indeed, the [ubiquitous|profound] interconnectedness of gene regulation networks has lead to a small revolution in our understanding of disease and complex traits [@Boyle2017-re].

Here, we use an alternative method of clustering to show that assortativity need not be the main driver of biologically meaningful clustering of complex traits. Using the Stochastic Block Model, which clusters genes so as to capture as much information on the network of interactions as possible, we find a large number of biologically relevant groups that could not be uncovered by methods that assume assortative modules. The blocks related to protein translation illustrate this nicely, with clear biological interpretation and practically no assortative modularity. This shifts our view of the structure of the relations between traits: instead of assuming the network is modular and clustering genes based on this assumption, we uncover clusters based on their information content about the gene network and ask if the resulting groups are assortative. Surprisingly, the answer is not always. We find assortative and non-assortative modules, and a marked difference in the overall modularity in the co-expression networks of two tissues in the same population.

There are also several practical advantages to the SBM. The non-parametric nature of the method means there are no free parameters to be optimized using heuristics, and other implementation choices, like the precise choice of edge weights, can be made in a principled fashion using the description length [@Peixoto2017-zw]. Also, the SBM finds a much larger number of modules that are guaranteed to be statistically supported, greatly improving the resolution of the clustering and allowing for more precise biological interpretation of the resulting blocks. The hierarchical nature of the model also allows for more or less coarse graining of the clustering, as illustrated by the signaling related block in fig. @fig:go_map. This very assortative block neatly decomposes into the elements involved in photoperception and neuronal signaling. 


- SBM recovered many more blocks and with more specific enrichment.
- No relation between assortativity and enrichment, several of the uncovered blocks are non-assortative.
- WGCNA and MMC fail to recover several clear blocks in SBM, like the neuro blocks in the head and the translation blocks in both tissues.
- SBM also reveals other biological insights, like the association between translation and cell respiration in the body.
- Imposing more structure than is warranted may be preventing us from extracting meaningful information from gene expression networks.
- Morphological traits being organized into modules can be interpreted as a reflexion of the very concrete structural and developmental constraints that lead to the formation and functioning of these individual body elements. No such clear imposition exists on gene expression, and the organization of genes into groups can happen through much more dynamic and varied mechanisms.

\newpage

# References





