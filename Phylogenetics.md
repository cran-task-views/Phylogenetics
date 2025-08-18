---
name: Phylogenetics
topic: Phylogenetics
maintainer: William Gearty, Brian O'Meara, Jacob Berv, Gustavo A. Ballen, Diniz Ferreira, Hilmar Lapp, Lars Schmitz, Martin R. Smith, Nathan S. Upham, Jonathan A. Nations
email: willgearty@gmail.com
version: 2025-02-17
source: https://github.com/cran-task-views/Phylogenetics/
---

## Overview

The history of life unfolds within a phylogenetic context, and phylogenetic trees (often shortened to "trees") are developed to represent this evolutionary history. Comparative phylogenetic methods are statistical approaches for analyzing historical patterns along such phylogenetic trees. This task view describes R packages that (i) facilitate the handling, manipulation and analysis of phylogenetic trees; (ii) implement comparative phylogenetic methods; (iii) apply phylogenetic methods to specific disciplines. This is an active research area and much of the information is subject to change. Many important packages are not on CRAN: either they were formerly on CRAN and were later archived (for example, if they failed to incorporate necessary changes as R is updated) or they are developed elsewhere and are not yet available on CRAN. Such packages may be found on GitHub, R-Forge, [Bioconductor](https://bioconductor.org/packages/release/BiocViews.html#___Phylogenetics), or authors' websites. At least ten packages start as phy\* in this domain, including two pairs of similarly named packages (phytools and phylotools, phylobase and phybase); users are encouraged to read and distinguish carefully between package names.

If you have any questions, feel free to reach out to the task view maintainers or the maintainers of specific packages. Questions may also be directed to the [R-SIG-Phylo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Phylo/) mailing-list after subscription.


## Scope

### Core packages

- `r pkg("ape", priority = "core")` implements the S3 phylo class which is commonly used to store phylogenetic trees in R. It is commonly used for reading, writing, and visualizing trees in the Newick/Phylip and NEXUS formats. It also has many functions for manipulating trees (e.g., rooting trees, dropping tips, randomly resolving polytomies), inferring trees (e.g., neighbour joining, bio-nj, and fast ME methods), and performing phylogenetic comparative analyses (e.g., reconstructing discrete and continuous characters, fitting basic models of trait evolution and diversification). It can also be used to generate random trees, pull in data from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/), and create lineage through time and correlogram plots.

- `r pkg("phylobase", priority = "core")` implements the S4 phylo4 class which combines phylogenetic trees and comparative data. While not used as commonly as the S3 phylo class, this new class is gaining traction among newer packages that implement phylogenetic comparative methods (e.g., `r pkg("adephylo")` and `r pkg("phylosignal")`).

- `r pkg("geiger", priority = "core")` implements a large suite of model fitting approaches for analyses of trait evolution and diversification. It is most commonly used to fit and compare various models of discrete and continuous trait evolution (e.g., Brownian motion, Ornstein-Uhlenbeck, Pagel's transforms, and models with trends). It also is commonly used to simulate phylogenies and the evolution of discrete and continuous characters. It also has several auxiliary functions that are often used by other packages.

- `r pkg("phytools", priority = "core")` has a constantly increasing range of functions for performing phylogenetic comparative analyses and visualizing (e.g., projecting into a morphospace), manipulating (e.g., branch length scaling and transformations, adding tips, finding subtrees), reading or writing, and even inferring phylogenetic trees and comparative data.

### Tasks

Packages within the task view fall within one or more of the following task categories:

1.  **Working with trees in R:** packages dedicated to the handling, manipulation, and visualization of phylogenetic data
2.  **Building trees in R:** packages for phylogenetic inference and tree simulation
3.  **Comparative phylogenetic methods:** packages for performing various comparative phylogenetic methods, including those dealing with trait evolution and diversification
4.  **Phylogenetics in other fields:** packages designed to perform field-specific phylogenetic analyses, including paleontology, community ecology, biogeography, and genetics
4.  **Other useful packages and miscellany:** packages that are useful for performing phylogenetic analyses, such as taxonomic matching


## Working with trees in R

### Getting trees into R

- `r pkg("phylobase")` and its lighter weight sibling `r pkg("rncl")` can use the [Nexus Class Library](http://ncl.sourceforge.net/) to read NEXUS, Newick, and other tree formats.
- `r pkg("treebase")` can search for and load trees from the online tree repository [TreeBASE](https://www.treebase.org/treebase-web/home.html).
- `r pkg("RNeXML")` can read, write, and process metadata for the [NeXML](http://www.nexml.org) format.
- `r pkg("TreeTools")` can read trees from external files in [TNT](https://cladistics.org/tnt/) format and NEXUS format, including extensions to the Nexus format not supported by `r pkg("ape")`, and metadata from [MorphoBank](https://morphobank.org/).
- `r pkg("ips")` can load trees from [BEAST](https://beast.community/), [MrBayes](http://nbisweden.github.io/MrBayes/), and other phylogenetics programs. This package can be used to parse the node support and other values from BEAST or MrBayes output.
- `r pkg("phylotate")` can read and write `r pkg("ape")`-compatible phylogenetic trees in NEXUS and Newick formats, while preserving annotations.
- `r pkg("phyext2")` can read and write various tree formats, including simmap formats.
- `r pkg("rotl")` can pull in a synthetic tree and individual study trees from the [Open Tree of Life](https://tree.opentreeoflife.org/) project.
- The `r bioc("treeio")` package can read trees in Newick, Nexus, New Hampshire eXtended format (NHX), jplace and Phylip formats and data output from BEAST, EPA, HyPhy, MrBayes, PAML, PHYLDOG, pplacer, r8s, RAxML and RevBayes.
- `r pkg("phylogram")` can convert Newick files into dendrogram objects.
- `r pkg("dendextend")` can manipulate such dendrogram objects.
- `r pkg("phytools")` can read and write trees in simple Newick and Nexus format, as well as `"simmap"` trees with an encoded discrete character.

### Tree manipulation

- `r pkg("phylobase")` has functions for traversing a tree (e.g., getting all descendants from a particular node specified by just two of its descendants).
- `r pkg("geiger")` can prune trees and data to an overlapping set of taxa. It can be also used to perform branch length scaling using ACDC; Pagel's (1999) lambda, delta and kappa parameters; and the Ornstein-Uhlenbeck alpha parameter (for ultrametric trees only). It can also be used to prune extinct taxa.
- `r pkg("TreeTools")` has functions to quantify and manipulate tree shape and balance, including the application of constraints; and to measure the phylogenetic information content of trees.
- `r pkg("Rogue")` identifies wildcard taxa, generating more informative summary trees.
- `r pkg("tidytree")` can convert a tree object to a tidy data frame and has other tidy approaches to manipulate tree data.
- `r pkg("evobiR")` can do fuzzy matching of names (to allow some differences).
- `r pkg("SigTree")` finds branches that are responsive to some treatment, while allowing correction for multiple comparisons.
- `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more.
- `r pkg("apex")` can handle multiple gene DNA alignments making their use and analysis for tree inference easier in `r pkg("ape")` and `r pkg("phangorn")`.
- `r pkg("aphid")` can weight sequences based on a phylogeny and can use hidden Markov models (HMMs) for a variety of purposes including multiple sequence alignment.
- `r pkg("phangorn")` and `r pkg("TreeSearch")` can perform tree rearrangements (NNI, SPR, and TBR).
- `r pkg("paleotree")` has functions for manipulating trees based on sampling issues that arise with fossil taxa as well as more universal transformations.
- `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more.
- `r pkg("castor")` can be used to manipulate extremely large trees (up to millions of tips).
- `r pkg("phytools")` can slice a tree at a pre-specified point, add taxa randomly to a tree, add species to genera, bind a single tip to a tree or two trees together, collapse clades on a tree using a clickable interface, perform midpoint rooting, paint a user-specified discrete character regime onto a tree to create a `"simmap"` object by various methods, convert a tree with a mapped character into a simple `"phylo"` object with unbranching nodes or a root edge into a single unbranching node, and other things.

### Tree visualization

- `r pkg("ape")`, `r pkg("adephylo")`, `r pkg("phylobase")`, `r pkg("phytools")`, `r pkg("ouch")`, and `r pkg("dendextend")` have functions for plotting trees; several of these have options for branch or taxon coloring based on some criterion (ancestral state, tree structure, etc.). In addition, `r pkg("phytools")` has substantial functionality to plot comparative data at the tips of the tree, graph the results of comparative analyses, and plot co-phylogenies.
- `r rforge("paleoPhylo")` and `r pkg("paleotree")` are specialized for drawing paleobiological phylogenies.
- `r github("heibl/viper")` can be used to annotate phylogenies with branch support, HPD intervals, and more.
- The popular R visualization package `r pkg("ggplot2")` can be extended by `r bioc("ggtree")` and `r bioc("ggtreeExtra")` to visualize phylogenies, and a geological timescale can be added using `r pkg("deeptime")`.
- `r pkg("strap")` can be used to add a geological timescale to a phylogeny, along with stratigraphic ranges.
- `r pkg("idendr0")` can be used to interactively explore trees (as dendrograms).
- `r pkg("phylocanvas")` is a widget for "htmlwidgets" that enables embedding of phylogenetic trees using the phylocanvas javascript library.
- `r pkg("ggmuller")` allows plotting a phylogeny along with frequency dynamics.
- `r pkg("RPANDA")` can be used to plot the spectral density and eigenvalues of a phylogeny.
- `r pkg("diversitree")` has an unexported function called "plot2.phylo()" which allows for the production of very lightweight PDF outputs of speciose trees (can be called via `diversitree:::plot2.phylo()`).

### Tree comparison

- `r pkg("distory")`, `r pkg("TreeDist")`, `r pkg("Quartet")`, `r pkg("TBRDist")`, `r pkg("phangorn")`, and `r pkg("phytools")` can compute distances between trees.
- `r pkg("TreeDist")` and `r pkg("treespace")` can plot and evaluate low-dimensional mappings of tree sets ("tree spaces").
- `r pkg("ape")` can compute tree-tree distances and also create a plot showing two trees with links between associated tips.
- `r pkg("dendextend")` can evaluate multiple measures comparing dendrograms.

### Phylogenetic summary statistics

- `r pkg("treestats")` can be used to calculate a wide collection of tree statistics, optimized for fast calculation.
- `r pkg("nLTT")` is specialized in calculating and visualising the nLTT statistic.
- `r pkg("castor")` contains fast calculation of the Gamma, Colless and Sackin statistic.
- `r pkg("phyloTop")` can be used to calculate a collection of pattern based summary statistics (e.g. ladders, stairs, cherries, pitchforks).
- `r pkg("treebalance")` can be used to calculate a collection of summary statistics focusing on measuring (im)balance.
- `r pkg("poweRbal")` provides an easy evaluation and comparison of tree shape statistics by estimating their power to differentiate between different tree models.
- `r pkg("RPANDA")` can compute Laplacian spectrum associated statistics.
- `r pkg("picante")` can compute community level summary statistics, such as mpd, mntd and psv.
- `r github("Leonardini/treeCentrality")` can compute several statistics inspired from network science.


## Tree building in R

### Phylogenetic inference

- `r pkg("phangorn")` can estimate trees using distance (e.g. UPGMA), parsimony, and likelihood.
- `r pkg("TreeSearch")` can identify most-parsimonious trees under parsimony, using the Brazeau et al. (2019) correction for inapplicable data, and includes a graphical user interface for detailed analysis of results.
- `r pkg("phyclust")` can cluster sequences.
- `r pkg("phytools")` can build trees using MRP supertree estimation and least squares.
- `r pkg("phylotools")` can build supermatrices for analyses in other software.
- `r pkg("EvoPhylo")` can be used to perform automated morphological character partitioning for bayesian phylogenetic analyses that are performed with [MrBayes](http://nbisweden.github.io/MrBayes/) and [BEAST2](https://www.beast2.org/). It can also be used to analyze the macroevolutionary parameter outputs from such analyses.
- `r bioc("fastreeR")` can be used to calculate distances, build phylogenetic trees, or perform hierarchical clustering between the samples of a VCF or FASTA file.
- `r github("uyedaj/rphenoscate")` facilitates analysis of "inapplicable" characters – such as morphological characters that are dependent on a parent trait – by constructing formal character hierarchies and the corresponding rate matrices (Tarasov 2023).

### Divergence times

- `r pkg("ape")` implements non-parametric rate smoothing (NPRS) and penalized likelihood.
- `r pkg("geiger")` can do congruification to stretch a source tree to match a specified standard tree.
- `r github("evolucionario/cladedate")` generates empirical calibration information from the fossil record.
- `r pkg("treedater")` implements various clock models, ways to assess confidence, and detecting outliers.
- `r pkg("phangorn")` can infer ultrametric and tipdated phylogenies with a strict clock model direct from sequences.
- `r github("dosreislab/bppr")` calibrates phylogenies from the program [BPP](https://github.com/bpp/bpp). A tutorial is available at [https://dosreislab.github.io/2018/08/31/bppr.html](https://dosreislab.github.io/2018/08/31/bppr.html)
- `r github("dosreislab/mcmc3r")` calculates the marginal likelihood in divergence time estimation using MCMCtree from the suite [PAML](https://github.com/abacus-gene/paml). It also calculates the block bootstrap for error estimation in marginal likelihood calculation.
- `r pkg("tbea")` allows to carry out multiple pre- (e.g. fitting densities to calibration quantiles, matrix and tree format conversion and matrix concatenation) and post-processing (e.g., summary of posterior tree samples, plotting prior vs. posterior estimates, measuring distribution similarity) tasks in Bayesian divergence time estimation. It has tools for summarizing collections of distributions (e.g. multiple estimates for the time of a biogeographic event, the origination time of a group where multiple estimates are available) which can be useful on their own or as a way to specify secondary caliibrations.

### Tree simulations

- `r pkg("TreeSim")` can be used to simulate trees using constant-rate birth-death with various constraints.
- `r pkg("phytools")` can simulate birth-death trees with various constraints, in both continuous and discrete time.
- `r pkg("geiger")` can be used to simulate trees under a birth-death process.
- `r pkg("DDD")` can be used to simulate trees under a diversity-dependent diversification process.
- `r pkg("paleotree")` can simulate fossil deposition, sampling, and the tree arising from this as well as trees conditioned on observed fossil taxa.
- `r pkg("FossilSim")` can be used to simulate fossil data on existing phylogenetic trees under mechanistic models of preservation and sampling.
- `r pkg("treats")` can be used to simulate birth-death phylogenetic trees and species traits jointly.
- `r pkg("TESS")` can simulate trees with time-dependent speciation and/or extinction rates, including mass extinctions.
- `r pkg("paleobuddy")` presents a flexible interface to simulate a wide array of user-defined diversification dynamics, including environmental-dependence.
- `r pkg("poweRbal")` provides a multitude of tree models to generate rooted binary trees with a given number of leaves.
- `r pkg("TreeSimGM")` can be used to simulate phylogenetic trees under general Bellman–Harris models with lineage-specific shifts of speciation and extinction, including simulating clade-dependent diversification processes.
- `r github("dosreislab/simclock")` simulates trees with branch lengths in number of substitutions per site under the relaxed clock models geometric Brownian motion (correlated rates) as well as independent lognormal rates.


## Phylogenetic comparative methods

### Ancestral state reconstruction

- `r pkg("ouch")` can be used to reconstruct root ancestral character states under Brownian motion or Ornstein-Uhlenbeck models, though ancestral states at the internal nodes are not.
- `r pkg("markophylo")` can fit a broad set of discrete character types with models that can incorporate constrained substitution rates, rate partitioning across sites, branch-specific rates, sampling bias, and non-stationary root probabilities.
- `r pkg("phytools")` can do ancestral character estimation for continuous and discrete characters under multiple models, including the threshold model from evolutionary quantitative genetics, as well as stochastic character mapping.
- `r pkg("Rphylopars")` can perform ancestral state reconstruction for datasets with multiple observations per species and/or missing data.
- `r pkg("TreeSearch")` can perform mapping of characters under parsimony, with an allowance for inapplicable data.
- `r pkg("castor")` can be used reconstruct continuous or discrete characters on extremely large trees.

### Trait evolution

- `r pkg("ape")`, `r pkg("picante")`, or `r pkg("caper")` can be used to calculate independent contrasts for continuous characters. `r pkg("caper")` also implements the brunch and crunch algorithms.
- `r pkg("geiger")` can be used to perform analyses of discrete trait evolution, including models of unequal rates or rates changing at a given instant of time, as well as Pagel's transformations.
- `r pkg("corHMM")` can be used to fit hidden Markov models of discrete character evolution which allow different transition rate classes on different portions of a phylogeny.
- `r pkg("phytools")` can be used to fit multiple models for both discrete and continuous character evolution. For instance, `r pkg("phytools")` fits a Brownian model with and without rate heterogeneity specified as regimes fixed by the user or estimated from the data itself. `r pkg("phytools")` can also be used to fit a range of discrete character evolution models, such as the the extended M*k* model, a heterogenous M*k* model with regime shifts, a polymorphic trait evolution model, a hidden-rates model, the threshold model, and others.
- `r pkg("geiger")`, `r pkg("paleotree")`, and `r pkg("motmot")` can be used to fit Brownian motion models.
- `r github("cran/RBrownie")` can fit multiple-rate Brownian motion models.
- `r pkg("geiger")` and `r pkg("OUwie")` can be used to investigate deviations from Brownian motion.
- `r pkg("mvMORPH")` can fit Brownian motion, early burst, ACDC, OU, and shift models to univariate or multivariate data.
- `r github("gilles-didier/cauphy")` models trait distribution using the Cauchy Process.
- `r pkg("geiger")`, `r pkg("motmot")`, `r pkg("ouch")`, `r pkg("slouch")`, `r pkg("surface")`, and `r pkg("OUwie")` can be used to fit Ornstein-Uhlenbeck (OU) models. `r pkg("ouch")` and `r pkg("slouch")` can implement models with multiple means, `r pkg("surface")` can implement models with multiple means using stepwise AIC, `r pkg("slouch")` can implement models with continuous covariates, and `r pkg("OUwie")` can implement models with multiple means, rates, and attraction values. Also see `r github("mongiardino/extendedSurface")` which combines the functionality of `r pkg("surface")` and `r pkg("OUwie")`.
- `r pkg("motmot")` can be used to fit continuous models that change rate or mode at specific time(s).
- `r pkg("Rphylopars")` can be used to fit continuous models such as those described above to datasets with multiple observations per species and/or missing data.
- `r pkg("geiger")` implements ANOVA's and MANOVA's in a phylogenetic context.
- `r pkg("ape")`, `r pkg("PHYLOGR")`, `r pkg("caper")`, and `r pkg("motmot")` implement traditional GLS methods (sensu Grafen or Martins).
- `r pkg("ape")` can be used to calculate phylogenetic autoregression (sensu Cheverud et al).
- `r pkg("ape")` and `r pkg("adephylo")` can be used to calculate phylogenetic autocorrelation (Moran's I).
- `r pkg("MCMCglmm")` can be used to assess correlation between traits using a GLMM.
- `r pkg("phylolm")` can fit phylogenetic linear regression and phylogenetic logistic regression models using a fast algorithm, making it suitable for large trees.
- `r pkg("brms")` can examine correlations between continuous and discrete traits, and can incorporate multiple measurements per species.
- `r pkg("metafor")` can perform meta-analyses accounting for phylogenetic structure.
- `r pkg("pmc")` evaluates the model adequacy of several trait models (from `r pkg("geiger")` and `r pkg("ouch")`) using Monte Carlo approaches.
- `r pkg("phyreg")` implements the Grafen (1989) phylogenetic regression.
- `r pkg("geomorph")` can do geometric morphometric analysis in a phylogenetic context.
- `r pkg("dispRity")` can be used to calculate disparity through time and perform other disparity-related analyses.
- `r pkg("MPSEM")` can predict features of one species based on information from related species using phylogenetic eigenvector maps.
- `r pkg("convevol")` and `r pkg("windex")` can both test for convergent evolution on a phylogeny.
- `r pkg("Claddis")` can be used to measure morphological diversity from discrete character data and evolutionary tempo on phylogenetic trees.
- `r pkg("do3PCA")` can be used to estimate probabilistic phylogenetic Principal Component Analysis (PCA), including methods to implement alternative models of trait evolution including Brownian motion (BM), Ornstein-Uhlenbeck (OU), Early Burst (EB), and Pagel's lambda.

### Trait simulations

- `r pkg("ouch")`, `r pkg("geiger")`, `r pkg("ape")`, `r pkg("picante")`, `r pkg("OUwie")`, `r pkg("caper")`, and `r pkg("phytools")` can be used to simulate continuous traits using Brownian motion.
- `r pkg("ouch")` and `r pkg("OUwie")` can be used to simulate continuous traits using the Hansen model (a form of the OU).
- `r pkg("geiger")` can be used to simulate continuous traits using a speciational model and discrete traits can be simulated using a continuous time Markov model (including models where rates change through time).
- `r pkg("phangorn")` can simulate DNA or amino acids.
- `r pkg("phytools")` can simulate discrete character evolution under multiple models.
- `r pkg("phylolm")` can simulate continuous or binary traits along a tree.
- `r pkg("Rphylopars")` can simulate data with missing observations.
- `r pkg("secsse")` can be used to simulate diversification models with a multistate observed trait and a hidden trait.
- `r pkg("treats")` can be used to simulate species traits and birth-death phylogenetic trees jointly.


### Diversification analysis

- `r pkg("ape")` and `r pkg("phytools")` can fit a simple birth-death model for when you have extant species only (sensu Nee et al. 1994), survival models, and goodness-of-fit tests (as applied to testing of models of diversification).
- `r pkg("TESS")` can calculate the likelihood of a tree under a model with time-dependent diversification, including mass extinctions.
- `r pkg("geiger")` can calculate net rates of diversification (sensu Magellon and Sanderson).
- `r pkg("diversitree")` implements the BiSSE method (Maddison et al. 1997) and later improvements (FitzJohn et al. 2009).
- `r pkg("hisse")` implements various hidden state diversification models, including HiSSE (Beaulieu and O'Meara 2016), GeoHiSSE (Caetano et al. 2018), MuHiSSE (Nakov et al. 2019), and MiSSE (trait-independent).
- `r pkg("caper")` can do the macrocaic test to evaluate the effect of a a trait on diversity.
- `r pkg("DDD")` implements maximum likelihood methods based on the diversity-dependent birth-death process to test whether speciation or extinction are diversity-dependent, as well as identifies key innovations and simulate a density-dependent process.
- `r pkg("PBD")` can calculate the likelihood of a tree under a protracted speciation model.
- `r pkg("phyloTop")` has functions for investigating tree shape, with special functions and datasets relating to trees of infectious diseases.
- `r pkg("RPANDA")` can be used to fit various diversification models to phylogenies, including time-dependent and environmental-dependent models.
- `r pkg("picante")` can be used to calculate various evolutionary distinctiveness measures, including the "equal splits" (ES) measure.
- `r pkg("castor")` can be used to estimate identifiable diversification rate parameters from trees (e.g., pulled rates of speciation).
- `r pkg("secsse")` can be used to fit diversification models with a multistate observed trait and a hidden trait.
- `r pkg("phytools")` can compute and visualize a lineages-through-time (LTT) plot, and calculate Pybus and Harvey's (2000) gamma statistic.
- `r pkg("CRABS")` features tools for exploring congruent phylogenetic birth-death models (see Louca and Pennell 2020).


## Phylogenetics in specific fields

### Morphometrics

- `r pkg("geomorph")` and `r pkg("RRPP")` may be used to evaluate evolutionary trends in multivariate phenotypes. Available methods include phylogenetic linear models (phylogenetic anova/regression, etc.), phylogenetic partial least squares, comparing rates of phenotypic evolution, phylogenetic integration, and phylogenetic modularity. Additionally, ordination approaches include both phylogenetic PCA, and phylogenetically-aligned components analysis (PACA).
- `r github("dosreislab/mcmc3r")` removes among-trait correlation using matrix shrinkage for estimation of branch lengths under Bronwnian motion and finally divergence time estimation in MCMCtree from the suite [PAML](https://github.com/abacus-gene/paml). A detailed tutorial can be found at [https://github.com/sabifo4/morpho](https://github.com/sabifo4/morpho).

### Time series and paleontology

- `r pkg("paleoTS")` can be used to analyze paleontological time series data using a likelihood-based framework for fitting and comparing models (using a model testing approach) of phyletic evolution (based on the random walk or stasis model).
- `r pkg("strap")` can do stratigraphic analysis of phylogenetic trees.
- `r github("rachelwarnock/fbdR")` can be used to estimate diversification rates from phylogenetic trees and fossil occurrence data.
- `r pkg("tbea")` has tools for estimating confidence intervals on stratigraphic end points and code for summarizing multiple distributions describing the same parameter.
- R offers a wealth of other options for general-purpose time series modeling, many of which are listed in the `r view("TimeSeries")` task view.

### Community and microbial ecology

- `r pkg("picante")`, `r pkg("vegan")`, `r pkg("SYNCSA")`, `r pkg("phylotools")`, `r pkg("caper")`, `r pkg("DAMOCLES")`, and `r pkg("phyloregion")` integrate several tools for using phylogenetics with community ecology.
- `r pkg("HMPTrees")` and `r pkg("GUniFrac")` provide tools for comparing microbial communities.
- `r pkg("betapart")` allows computing pair-wise dissimilarities (distance matrices) and multiple-site dissimilarities, separating the turnover and nestedness-resultant components of taxonomic (incidence and abundance based), functional and phylogenetic beta diversity.
- `r pkg("phyloregion")` extends  `r pkg("betapart")` to allow sparse community matrices allowing larger datasets.
- `r pkg("adiv")` can calculate various indices of biodiversity including species, functional and phylogenetic diversity, as well as alpha, beta, and gamma diversities.
- `r pkg("entropart")` can measure and partition diversity based on Tsallis entropy as well as calculate alpha, beta, and gamma diversities.
- `r pkg("metacoder")` provides functions for handling large taxonomic data sets, like those generated from modern high-throughput sequencing, like metabarcoding.
- `r bioc("phyloseq")` provides a set of classes and tools to facilitate the import, storage, analysis, and graphical display of microbiome census data.

### Phyloclimatic modeling

- `r github("heibl/phyloclim")` integrates several tools for phyloclimatic modeling.

### Phylogeography and biogeography

- `r pkg("diversitree")` implements the GeoSSE method for diversification analyses based on two areas.
- `r github("GabrielNakamura/Herodotools")` can be used to perform a wide variety of biogeographical macroevolutionary analyses, including exploring spatial biodiversity patterns, conducting ancestral area reconstruction, and performing evoregion classification.
- `r pkg("epm")` (EcoPhyloMapper) can be used to calculate various morphological and phylogenetic community metrics across geography.
- `r github("nmatzke/BioGeoBEARS")` allows probabilistic inference of both historical biogeography (ancestral geographic ranges on a phylogeny) as well as comparison of different models of range evolution.

### Epidemiology

See the `r view("Epidemiology")` task view for details about packages useful for epidemiology, including phylogenetic epidemiology.

### Omics

- `r pkg("aphylo")` implements a parsimonious evolutionary model to analyze and predict gene-functional annotations in phylogenetic trees.
- `r pkg("CALANGO")` can be used to search for annotation terms (e.g., Pfam IDs, GO terms or superfamilies) associated with a quantitative/rank variable.
- `r github("hr1912/TreeExp")` can be used to perform comparative analyses of gene expression in a phylogenetic context.
- See the `r view("Omics")` task view for details about other useful packages.

### Gene tree--species tree and species delimitation

- `r rforge("splits")` uses a gene tree to infer species limits based on GMYC (Generalized Mixed Yule Coalescent).
- `r github("emanuelmfonseca/P2C2M.GMYC")` can identify model violations under a GMYC model.
- `r pkg("MSCquartets")` provides methods for analyzing and using quartets displayed on a collection of gene trees, primarily to make inferences about the species tree or network under the multi-species coalescent model.
- `r github("dosreislab/bppr")` can prepare the control files for doing model selection for comparing competing species trees using [BPP](https://github.com/bpp/bpp). A tutorial is available at [https://dosreislab.github.io/2018/08/31/bppr.html](https://dosreislab.github.io/2018/08/31/bppr.html).


## Other useful packages and miscellany

### Taxonomy

- `r pkg("taxize")` can interact with a suite of web APIs for taxonomic tasks, such as verifying species names, getting taxonomic hierarchies, and verifying name spelling.
- `r pkg("evobiR")` contains functions for making a tree at higher taxonomic levels, downloading a taxonomy tree from NCBI or ITIS, and various other miscellaneous functions (simulations of character evolution, calculating D-statistics, etc.).

### Interactions with other programs

- `r pkg("ape")` can call [PhyML](http://www.atgc-montpellier.fr/phyml/), [Clustal](http://www.clustal.org/), [T-Coffee](https://tcoffee.crg.eu/), and [Muscle](https://www.drive5.com/muscle/) through various functions.
- `r pkg("geiger")` can call PATHd8 through its congruify function.
- `r pkg("ips")` wraps several phylogenetic software for sequence alignment, masking of sequence alignments, and estimation of phylogenies and ancestral character states, including MrBayes, Beast, [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/), [PartitionFinder](https://www.robertlanfear.com/partitionfinder/), and [MAFFT](https://mafft.cbrc.jp/alignment/software/), allowing their easy use from within R.
- `r pkg("beastier")` can call [BEAST2](https://www.beast2.org/) to perform phylogenetic analyses, `r pkg("beautier")` can generate XML input files for BEAST2 (like [BEAUti](https://www.beast2.org/beauti/)), and `r pkg("tracerer")` can be used to parse and analyze BEAST2 output files (like [Tracer](https://github.com/beast-dev/tracer/)). `r pkg("babette")` is a wrapper for all of these packages. `r pkg("mcbette")` allows to do Bayesian model comparison over some site and clock models using `r pkg("babette")`.
- `r github("liamrevell/Rphylip")` wraps [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) , a broad variety of programs for tree inference under parsimony, likelihood, and distance, bootstrapping, character evolution, and more.
- `r pkg("BoSSA")` can use information from various tools to place a query sequence into a reference tree.
- `r pkg("BAMMtools")` is an interface to the BAMM program to allow visualization of rate shifts, comparison of diversification models, and other functions.
- `r pkg("Revticulate")` can be used to interact with [RevBayes](https://revbayes.github.io/) from within R, while `r pkg("RevGadgets")` can be used to process the output generated by RevBayes.
- `r github("dosreislab/bppr")` can prepare the control files for doing model selection for comparing competing species trees using [BPP](https://github.com/bpp/bpp). A tutorial is available at [https://dosreislab.github.io/2018/08/31/bppr.html](https://dosreislab.github.io/2018/08/31/bppr.html)
- `r github("dosreislab/mcmc3r")` can prepare the control files for carrying out divergence time estimation using MCMCtree from the suite [PAML](https://github.com/abacus-gene/paml). It also generates morphological alignments in phylip format for using continuous trait models in divergence time estimation in MCMCtree.
- `r pkg("tbea")` has code for post-analysis summarization and plotting of trace files from Bayesian phylogenetic programs such as Beast2, MrBayes, RevBayes, and MCMCTree.


## References

-   Beaulieu, J.M. and O'Meara, B.C., 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology, 65(4): 583-601. `r doi("10.1093/sysbio/syw022")`.
-   Borregaard, M.K., Rahbek, C., Fjeldsaa, J., Parra, J.L., Whittaker, R.J. and Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution 5(11): 1225-1235. `r doi("10.1111/2041-210X.12283")`.
-   Brazeau, M.D., Guillerme, T. and Smith, M.R. 2019. An algorithm for morphological phylogenetic analysis with inapplicable data. Systematic Biology, 68:619--631. `r doi("10.1093/sysbio/syy083")`.
-   Butler M.A., King A.A. 2004 Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164, 683-695. `r doi("10.1086/426002")`.
-   Caetano, D.S., B.C. O'Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, 72:2308-2324. `r doi("10.1111/evo.13602")`.
-   Cheverud J.M., Dow M.M., Leutenegger W. 1985 The quantitative assessment of phylogenetic constraints in comparative analyses: Sexual dimorphism in body weight among primates. Evolution 39, 1335-1351. `r doi("10.1111/j.1558-5646.1985.tb05699.x")`.
-   FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58: 595-611. `r doi("10.1093/sysbio/syp067")`.
-   Garland T., Harvey P.H., Ives A.R. 1992 Procedures for the analysis of comparative data using phylogenetically independent contrasts. Systematic Biology 41, 18-32. `r doi("10.1093/sysbio/41.1.18")`.
-   Hansen T.F. 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution 51: 1341-1351. `r doi("10.1111/j.1558-5646.1997.tb01457.x")`.
-   Louca, S. and Pennell, M. W. 2020. Extant timetrees are consistent with a myriad of diversification histories. Nature, 580(7804), 502-505. `r doi("10.1038/s41586-020-2176-1")`.
-   Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary character's effect on speciation and extinction. Systematic Biology 56: 701--710. `r doi("10.1080/10635150701607033")`.
-   Magallon S., Sanderson, M.J. 2001. Absolute Diversification Rates in Angiosperm Clades. Evolution 55(9):1762-1780. `r doi("10.1111/j.0014-3820.2001.tb00826.x")`.
-   Moore, B.R., Chan, K.M.A., Donoghue, M.J. (2004) Detecting diversification rate variation in supertrees. In Bininda-Emonds ORP (ed) Phylogenetic Supertrees: Combining Information to Reveal the Tree of Life, Kluwer Academic pgs 487-533. `r doi("10.1007/978-1-4020-2330-9_23")`.
-   Nakov, T., Beaulieu, J.M., and Alverson, A.J. 2019. Diatoms diversify and turn over faster in freshwater than marine environments. Evolution, 73: 2497-2511. `r doi("10.1111/evo.13832")`.
-   Nee S., May R.M., Harvey P.H. 1994. The reconstructed evolutionary process. Philosophical Transactions of the Royal Society of London Series B Biological Sciences 344: 305-311. `r doi("10.1098/rstb.1994.0068")`.
-   Pagel M. 1999. Inferring the historical patterns of biological evolution. Nature 401, 877-884. `r doi("10.1038/44766")`.
-   Pybus O.G., Harvey P.H. 2000. Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London Series B Biological Sciences 267, 2267-2272. `r doi("10.1098/rspb.2000.1278")`.
-   Tarasov, S. 2023. New phylogenetic Markov models for inapplicable morphological characters. Systematic Biology, syad005. `r doi("10.1093/sysbio/syad005")`.
