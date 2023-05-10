Here we have used the Stochastic Block Model, a clustering method that does not rely on modularity maximization, to explore the organization of gene co-expression networks in female Drosophila melanogaster.
The SBM, in contrast with other methods explored here, was able to a) cluster all genes into blocks, b) identify blocks with both, high resolution (few genes per block) and high functional content (significant GO associations), and c) identify blocks that are modular (higher within- than between-block correlation) as well as non-modular.
The high-resolution and hierarchical manner in which the transcriptomes were resolved offers the complex trait community a powerful resource to infer the function of genes using a guilty-by-association approach, as it allows for the use of contextual information given by the hierarchy when making precise inferences about unannotated genes in small blocks.


Besides the better performance of the SBM in finding gene clusters, when compared to other community detection methods, there are also several practical advantages offered by the SBM.
For example, in WGCNA there are three clustering (?) parameters that need to be manually tuned, and the specific choices can drastically change the number of genes that are clustered and the number and size of modules.
In the SBM, on the other hand, the clustering procedure is completely parameter-free, and choices regarding how to model the weights between edges can be made by selecting the model with the shortest description length (Peixoto, 2017).
While there are some heuristics for tuning the WGCNA parameters, they rely heavily on the expectation that gene co-expression networks should be approximately scale-free (Bergmann et al., 2004; Dong & Horvath, 2007; Jeong et al., 2000), but, despite its popularity, this expectation might be unwarranted (Khanin & Wit, 2006; Stumpf et al., 2005).
Therefore, using the SBM offers a significant advantage in situations where we lack a priori information to tune such parameters.

Even with optimal parameters, WGCNA fails to assign a substantial proportion of genes to any module.
This is because, while WGCNA is efficient in finding hub genes, if some functional gene group does not have a hub or has low average similarities, this group will never be identified.
SBM, on the contrary, finds a much larger number of modules that are guaranteed to be statistically supported, greatly improving the resolution of the clustering and allowing for more precise biological interpretation of the resulting blocks.
And, in addition, the posterior distribution of block assignments can be used to quantify the confidence level of the grouping of particular genes or blocks.
While we have not explored the statistical side of the SBM here, there is ample space for using the SBM as a theoretically well-supported statistical tool in gene expression clustering.
Biology has a history of borrowing methods and insight from other fields like xxxx, xxx, and xxxx for answering biological questions (O.
Mason & Verwoerd, 2007; Proulx et al., 2005; Radde & Hütt, 2016).
Network analysis has come a long way since modularity maximization (Peixoto, 2021), and the gene expression and modularity community should continue this rich tradition by porting recent advancements from network theory, like the SBM, into our toolkit.

One aspect we did not explore here is the estimation of the gene co-expression network itself, before any attempt at finding communities.
Both types of methods used weighted networks: fully connected ones for the WGCNA and MMC pipelines (as per these methods’ suggested workflow) and a sparser network for the SBM model fitting, due to computational constraints.
Estimating these weights (gene expression correlations) is an error-prone process, as we are estimating many more weights than we have measured individuals, leading to potentially poor estimates (Schäfer & Strimmer, 2005).
While the procedures we used here are commonplace, there are more principled ways of building co-expression networks (see Peel et al., 2022 for a recent perspective), and this is an aspect of the usual transcriptomics workflow that could potentially see massive improvements in the near future.
Methods like the graphical lasso have been used in this context (Lingjærde et al., 2021; Lyu et al., 2018; Seal et al., 2023), and the expectation is that, when compared to fully connected or thresholded networks, these inferred networks should provide much better estimates of gene-gene connections and weights.
Additionally, it is possible to combine community detection via the SBM with network inference, simultaneously using information about community structure to inform the network inference and vice-versa (Peixoto, 2019).

Beyond the methodological and practical advantages discussed above, the fact that the SBM does not find gene clusters by attempting to maximize their modularity has major implications for our understanding of biological networks.
Here we find that D.
melanogaster transcriptomes are organized into modular as well as non-modular gene clusters.
The latter, however, could not have been identified by methods that assume assortative modules.
The possibility of quantifying, in a continuous scale, the degree of modularity (assortativity) of each gene block allowed us to compare the gene co-expression networks derived from head and body tissue, and uncover marked differences in their overall degree of modularity.
This opens the possibility of expanding this comparison to different cell types, organs, and even species to get a comprehensive understanding of how modular are indeed biological networks.


These results warrant a reassessment of the assumption that gene co-expression networks are modular.
Modularity, understood as the relative independence between parts of complex traits (e.g., organs, networks), is often invoked to explain the evolvability of complex phenotypes and has functioned as a unifying concept at several levels of organization with great success (Melo et al., 2016; Wagner et al., 2007; Zelditch & Goswami, 2021).
A modular organization allows organisms to maintain their function unchanged by coordinating simultaneous evolutionary changes in all related traits while keeping unrelated traits undisturbed (Ancel & Fontana, 2000; Cheverud, 1996; Wagner & Altenberg, 1996; Wagner & Zhang, 2011).
A large part of the literature on modularity developed in the context of morphological traits, where it is easy to understand the structural and developmental constraints that a modular organization helps to overcome, and therefore allowing change while maintaining a functionality (Marcucio et al., 2011; Shirai & Marroig, 2010).

However, no such clear structural and physical constraints exist on gene expression, and the interaction between groups of genes can happen through much more dynamic and varied mechanisms.
While we might expect related genes to be co-expressed and therefore highly correlated, non-linear phenomena can lead to a complete decoupling of the expression levels of co-expressed genes.
For example, the effect of gene A on gene B could have a saturation point after which increasing expression of gene A no longer leads to higher levels of gene B, and no correlation is detected in this regime, even if the genes are co-expressed.
The marked difference in the level of modularity across the two tissues in our samples illustrates just how variable modularity can be, even within the same species, sex, and population.
Furthermore, modularity is not a necessary feature of biological organization (even in the case of evolvability, see Pavlicev & Hansen, 2011; Roseman et al., 2009), and only searching for modularity can blind us to alternative organizations, as we have shown.
Indeed, the realization of the profound interconnectedness of gene regulation networks has led to a small revolution in our understanding of disease and complex traits (Boyle et al., 2017).
The very high dimensionality of gene co-expression networks also allows for genes to be similar in ways that do not lead to high correlations.
For example, two genes might be connected to the same genes in different modules, but not among themselves.
This similarity would likely be missed by modularity maximization because these genes would not form a classic assortative unit.
Meanwhile, the SBM would correctly identify these genes connecting two modules as being similar due to their shared connectivity pattern.


Here we find that non-modular blocks are widespread in gene co-expression networks, and that the evidence for their functional relevance is as strong as for modular blocks.
This highlights the need to incorporate other sources of information, beyond modularity, when exploring biological networks.
More studies using methods that don’t rely on modularity maximization will be needed to determine whether there are general patterns of non-modular organization.
For example, here we find that, despite the differences in gene clusters between body and head, the non-modular blocks tend to be associated with protein translation.
Will this emerge as a generality of transcriptomes?