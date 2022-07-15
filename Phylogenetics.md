---
name: Phylogenetics
topic: Phylogenetics packages in R
maintainer: Brian O'Meara, William Gearty
email: omeara.brian@gmail.com, willgearty@gmail.com
version: 2022-07-11
source: https://github.com/cran-task-views/Spatial/
---

## Overview

The history of life unfolds within a phylogenetic context. Comparative phylogenetic methods are statistical approaches for analyzing historical patterns along phylogenetic trees. This task view describes R packages that 1) are useful for the handling and manipulation of phylogenetic trees in R and/or 2) implement a variety of different comparative phylogenetic methods. This is an active research area and much of the information is subject to change. One thing to note is that many important packages are not on CRAN: either they were formerly on CRAN and were later archived (for example, if they failed to incorporate necessary changes as R is updated) or they are developed elsewhere and have not been put on CRAN yet. Such packages may be found on GitHub, R-Forge, or authors' websites.

Packages are grouped into the following categories:

1.  **Working with trees in R:** packages dedicated to the handling, manipulation, and visualization of phylogenetic data
2.  **Building trees in R:** packages for phylogenetic inference and tree simulation
3.  **Comparative phylogenetic methods:** packages for performing various comparative phylogenetic methods, including those dealing with trait evolution and diversification
4.  **Phylogenetics in other fields:** packages designed to perform field-specific phylogenetic analyses, including paleontology, community ecology, biogeography, and genetics
4.  **Other useful packages and miscellany:** packages that are useful for performing phylogenetic analyses, such as taxonomic matching

## Working with trees in R

### Getting trees into R

Trees in R are usually stored in the S3 phylo class (implemented in `r pkg("ape", priority = "core")`), though the S4 phylo4 class (implemented in `r pkg("phylobase")`) is also available. `r pkg("ape")` can read trees from external files in newick format (sometimes popularly known as phylip format) or NEXUS format. It can also read trees input by hand as a newick string (i.e., "(human,(chimp,bonobo));"). `r pkg("phylobase")` and its lighter weight sibling `r pkg("rncl")` can use the [Nexus Class Library](http://ncl.sourceforge.net/) to read NEXUS, Newick, and other tree formats. `r pkg("treebase")` can search for and load trees from the online tree repository TreeBASE, `r pkg("rdryad")` can pull data from the online data repository Dryad. `r pkg("RNeXML")` can read, write, and process metadata for the [NeXML](http://www.nexml.org) format. PHYLOCH can load trees from BEAST, MrBayes, and other phylogenetics programs (PHYLOCH is only available from the author's [website](http://www.christophheibl.de/Rpackages.html) ). `r pkg("phyext2")` can read and write various tree formats, including simmap formats. `r pkg("rotl")` can pull in a synthetic tree and individual study trees from the Open Tree of Life project. The `r bioc("treeio")` package can read trees in Newick, Nexus, New Hampshire eXtended format (NHX), jplace and Phylip formats and data output from BEAST, EPA, HyPhy, MrBayes, PAML, PHYLDOG, pplacer, r8s, RAxML and RevBayes. `r pkg("phylogram")` can convert Newick files into dendrogram objects. `r pkg("brranching")` can fetch phylogenies from online repositories, including [phylomatic](http://phylodiversity.net/phylomatic/) .

### Utility functions

These packages include functions for manipulating trees or associated data. `r pkg("ape")` has functions for randomly resolving polytomies, creating branch lengths, getting information about tree size or other properties, pulling in data from GenBank, and many more. `r pkg("phylobase")` has functions for traversing a tree (i.e., getting all descendants from a particular node specified by just two of its descendants). `r pkg("geiger")` can prune trees and data to an overlapping set of taxa. `r pkg("tidytree")` can convert a tree object in to a tidy data frame and has other tidy approaches to manipulate tree data. `r pkg("evobiR")` can do fuzzy matching of names (to allow some differences). `r pkg("SigTree")` finds branches that are responsive to some treatment, while allowing correction for multiple comparisons. `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more. `r pkg("apex")` can handle multiple gene DNA alignments making their use and analysis for tree inference easier in `r pkg("ape")` and `r pkg("phangorn")`. `r pkg("aphid")` can weight sequences based on a phylogeny and can use hidden Markov models (HMMs) for a variety of purposes including multiple sequence alignment.

### Tree manipulation

Branch length scaling using ACDC; Pagel's (1999) lambda, delta and kappa parameters; and the Ornstein-Uhlenbeck alpha parameter (for ultrametric trees only) are available in `r pkg("geiger")`. `r pkg("phytools")` also allows branch length scaling, as well as several tree transformations (adding tips, finding subtrees). Rooting, resolving polytomies, dropping of tips, setting of branch lengths including Grafen's method can all be done using `r pkg("ape")`. Extinct taxa can be pruned using `r pkg("geiger")`. `r pkg("phylobase")` offers numerous functions for querying and using trees (S4). Tree rearrangements (NNI and SPR) can be performed with `r pkg("phangorn")`. `r pkg("paleotree")` has functions for manipulating trees based on sampling issues that arise with fossil taxa as well as more universal transformations. `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more.

### Tree visualization

User trees can be plotted using `r pkg("ape")`, `r pkg("adephylo")`, `r pkg("phylobase")`, `r pkg("phytools")`, `r pkg("ouch")`, and `r pkg("dendextend")`; several of these have options for branch or taxon coloring based on some criterion (ancestral state, tree structure, etc.). `r rforge("paleoPhylo")` and `r pkg("paleotree")` are specialized for drawing paleobiological phylogenies. Trees can also be examined (zoomed) and viewed as correlograms using `r pkg("ape")`. Ancestral state reconstructions can be visualized along branches using `r pkg("ape")` and `r pkg("paleotree")`. `r pkg("phytools")` can project a tree into a morphospace. `r pkg("BAMMtools")` can visualize rate shifts calculated by BAMM on a tree. The popular R visualization package `r pkg("ggplot2")` can be extended by `r github("GuangchuangYu/ggtree")` to visualize phylogenies, and a geological timescale can be added using `r pkg("deeptime")`. `r pkg("strap")` can also be used to add a geological timescale to a phylogeny, along with stratigraphic ranges. Trees can also be interactively explored (as dendrograms) using `r pkg("idendr0")`. `r pkg("phylocanvas")` is a widget for "htmlwidgets" that enables embedding of phylogenetic trees using the phylocanvas javascript library. `r pkg("ggmuller")` allows plotting a phylogeny along with frequency dynamics. `r pkg("RPANDA")` can be used to plot the spectral density and eigenvalues of a phylogeny.

### Tree comparison

Tree-tree distances can be evaluated, and used in additional analyses, in `r pkg("distory")` and `r pkg("Rphylip")`. `r pkg("ape")` can compute tree-tree distances and also create a plot showing two trees with links between associated tips. `r pkg("kdetrees")` implements a non-parametric method for identifying potential outlying observations in a collection of phylogenetic trees, which could represent inference problems or processes such as horizontal gene transfer. `r pkg("dendextend")` can evaluate multiple measures comparing dendrograms.

## Tree building in R

### Phylogenetic inference

Neighbour joining, bio-nj and fast ME methods of phylogenetic reconstruction are all implemented in the package `r pkg("ape")`. `r pkg("phangorn")` can estimate trees using distance (e.g. UPGMA), parsimony, and likelihood. `r pkg("phyclust")` can cluster sequences. `r pkg("phytools")` can build trees using MRP supertree estimation and least squares. `r pkg("phylotools")` can build supermatrices for analyses in other software. `r pkg("pastis")` can use taxonomic information to make constraints for Bayesian tree searches. For more information on importing sequence data, see the `r view("Genetics")` task view; `r pkg("pegas")` may also be of use. `r pkg("EvoPhylo")` can be used to perform automated morphological character partitioning for bayesian phylogenetic analyses that are performed with [MrBayes](http://nbisweden.github.io/MrBayes/) and [BEAST2](https://www.beast2.org/). It can also be used to analyze the macroevolutionary parameter outputs from such analyses. `r pkg("Revticulate")` can be used to interact with [RevBayes](https://revbayes.github.io/) from within R, while `r pkg("RevGadget")` can be used to process the output generated by RevBayes.

### Divergence times

Non-parametric rate smoothing (NPRS) and penalized likelihood can be implemented in `r pkg("ape")`. `r pkg("geiger")` can do congruification to stretch a source tree to match a specified standard tree. `r pkg("treedater")` implements various clock models, ways to assess confidence, and detecting outliers. `r pkg("phangorn")` can infer ultrametric and tipdated phylogenies with a strict clock model direct from sequences.  

### Tree simulations

Trees can be simulated using constant-rate birth-death with various constraints in `r pkg("TreeSim")` and a birth-death process in `r pkg("geiger")`. Random trees can be generated in `r pkg("ape")` by random splitting of edges (for non-parametric trees) or random clustering of tips (for coalescent trees). `r pkg("paleotree")` can simulate fossil deposition, sampling, and the tree arising from this as well as trees conditioned on observed fossil taxa. `r pkg("FossilSim")` can be used to simulate fossil data on existing phylogenetic trees under mechanistic models of preservation and sampling. `r pkg("TESS")` can simulate trees with time-dependent speciation and/or extinction rates, including mass extinctions. `r pkg("paleobuddy")` presents a flexible interface to simulate a wide array of user-defined diversification dynamics, including environmental-dependence.

## Comparative phylogenetic methods

### Ancestral state reconstruction

Continuous characters can be reconstructed using maximum likelihood, generalised least squares or independent contrasts in `r pkg("ape")`. Root ancestral character states under Brownian motion or Ornstein-Uhlenbeck models can be reconstructed in `r pkg("ouch")`, though ancestral states at the internal nodes are not. Discrete characters can be reconstructed using a variety of Markovian models that parameterize the transition rates among states using `r pkg("ape")`. `r pkg("markophylo")` can fit a broad set of discrete character types with models that can incorporate constrained substitution rates, rate partitioning across sites, branch-specific rates, sampling bias, and non-stationary root probabilities. `r pkg("phytools")` can do stochastic character mapping of traits on trees. Ancestral state reconstruction for datasets with multiple observations per species and/or missing data can be performed using `r pkg("Rphylopars")`.

### Trait evolution

Independent contrasts for continuous characters can be calculated using `r pkg("ape")`, `r pkg("picante")`, or `r pkg("caper")` (which also implements the brunch and crunch algorithms). Analyses of discrete trait evolution, including models of unequal rates or rates changing at a given instant of time, as well as Pagel's transformations, can be performed in `r pkg("geiger")`. Brownian motion models can be fit in `r pkg("geiger")`, `r pkg("ape")`, and `r pkg("paleotree")`. Deviations from Brownian motion can be investigated in `r pkg("geiger")` and `r pkg("OUwie")`. `r pkg("mvMORPH")` can fit Brownian motion, early burst, ACDC, OU, and shift models to univariate or multivariate data. Ornstein-Uhlenbeck (OU) models can be fitted in `r pkg("geiger")`, `r pkg("ape")`, `r pkg("ouch")` (with multiple means), `r pkg("surface")` (with multiple means using stepwise AIC), and `r pkg("OUwie")` (with multiple means, rates, and attraction values). Also see `r github("mongiardino/extendedSurface")` for a combines the functionality of `r pkg("surface")` and `r pkg("OUwie")`. `r pkg("geiger")` fits only single-optimum models. Other continuous models, including Pagel's transforms and models with trends, can be fit with `r pkg("geiger")`. Continuous models such as those described above can be fit to datasets with multiple observations per species and/or missing data using `r pkg("Rphylopars")`. ANOVA's and MANOVA's in a phylogenetic context can also be implemented in `r pkg("geiger")`. Multiple-rate Brownian motion can be fit in `r github("cran/RBrownie")`. Traditional GLS methods (sensu Grafen or Martins) can be implemented in `r pkg("ape")`, `r pkg("PHYLOGR")`, or `r pkg("caper")`. Phylogenetic autoregression (sensu Cheverud et al) and Phylogenetic autocorrelation (Moran's I) can be implemented in `r pkg("ape")` or\--if you wish the significance test of Moran's I to be calculated via a randomization procedure\--in `r pkg("adephylo")`. Correlation between traits using a GLMM can also be investigated using `r pkg("MCMCglmm")`. `r pkg("phylolm")` can fit phylogenetic linear regression and phylogenetic logistic regression models using a fast algorithm, making it suitable for large trees. `r pkg("brms")` can examine correlations between continuous and discrete traits, and can incorporate multiple measurements per species. `r pkg("phytools")` can also investigate rates of trait evolution and do stochastic character mapping. `r pkg("metafor")` can perform meta-analyses accounting for phylogenetic structure. `r pkg("pmc")` evaluates the model adequacy of several trait models (from `r pkg("geiger")` and `r pkg("ouch")`) using Monte Carlo approaches. `r pkg("phyreg")` implements the Grafen (1989) phyglogenetic regression. `r pkg("geomorph")` can do geometric morphometric analysis in a phylogenetic context. Disparity through time, and other disparity-related analyses, can be performed with `r pkg("dispRity")`. `r pkg("MPSEM")` can predict features of one species based on information from related species using phylogenetic eigenvector maps. `r pkg("Rphylip")` wraps [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) which can do independent contrasts, the threshold model, and more. `r pkg("convevol")` and `r pkg("windex")` can both test for convergent evolution on a phylogeny. `r pkg("Claddis")` can be used to measure morphological diversity from discrete character data and evolutionary tempo on phylogenetic trees.

### Trait simulations

Continuous traits can be simulated using brownian motion in `r pkg("ouch")`, `r pkg("geiger")`, `r pkg("ape")`, `r pkg("picante")`, `r pkg("OUwie")`, and `r pkg("caper")`, the Hansen model (a form of the OU) in `r pkg("ouch")` and `r pkg("OUwie")` and a speciational model in `r pkg("geiger")`. Discrete traits can be simulated using a continuous time Markov model in `r pkg("geiger")`. `r pkg("phangorn")` can simulate DNA or amino acids. Both discrete and continuous traits can be simulated under models where rates change through time in `r pkg("geiger")`. `r pkg("phytools")` can simulate discrete characters using stochastic character mapping. `r pkg("phylolm")` can simulate continuous or binary traits along a tree. Data with missing observations can be simulated using `r pkg("Rphylopars")`.

### Diversification analysis

Lineage through time plots can be done in `r pkg("ape")`. A simple birth-death model for when you have extant species only (sensu Nee et al. 1994) can be fitted in ape as can survival models and goodness-of-fit tests (as applied to testing of models of diversification). `r pkg("TESS")` can calculate the likelihood of a tree under a model with time-dependent diversification, including mass extinctions. Net rates of diversification (sensu Magellon and Sanderson) can be calculated in `r pkg("geiger")`. `r pkg("diversitree")` implements the BiSSE method (Maddison et al. 1997) and later improvements (FitzJohn et al. 2009). `r pkg("TreePar")` estimates speciation and extinction rates with models where rates can change as a function of time (i.e., at mass extinction events) or as a function of the number of species. `r pkg("caper")` can do the macrocaic test to evaluate the effect of a a trait on diversity. `r pkg("apTreeshape")` also has tests for differential diversification (see [description](https://doi.org/10.1093/bioinformatics/bti798) ). `r pkg("iteRates")` can identify and visualize areas on a tree undergoing differential diversification. `r pkg("DDD")` can fit density dependent models as well as models with occasional escape from density-dependence. `r pkg("BAMMtools")` is an interface to the BAMM program to allow visualization of rate shifts, comparison of diversification models, and other functions. `r pkg("DDD")` implements maximum likelihood methods based on the diversity-dependent birth-death process to test whether speciation or extinction are diversity-dependent, as well as identifies key innovations and simulate a density-dependent process. `r pkg("PBD")` can calculate the likelihood of a tree under a protracted speciation model. `r pkg("phyloTop")` has functions for investigating tree shape, with special functions and datasets relating to trees of infectious diseases. `r pkg("RPANDA")` can be used to fit various diversification models to phylogenies, including time-dependent and environmental-dependent models.

## Phylogenetics in specific fields

### Time series and paleontology

Paleontological time series data can be analyzed using a likelihood-based framework for fitting and comparing models (using a model testing approach) of phyletic evolution (based on the random walk or stasis model) using `r pkg("paleoTS")`. `r pkg("strap")` can do stratigraphic analysis of phylogenetic trees. R offers a wealth of other options for general-purpose time series modeling, many of which are listed in the `r view("TimeSeries")` task view. `r github("rachelwarnock/fbdR")` can be used to estimate diversification rates from phylogenetic trees and fossil occurrence data.

### Community and microbial ecology

`r pkg("picante")`, `r pkg("vegan")`, `r pkg("SYNCSA")`, `r pkg("phylotools")`, `r pkg("PCPS")`, `r pkg("caper")`, `r pkg("DAMOCLES")`, `r pkg("phyloregion")` integrate several tools for using phylogenetics with community ecology. `r pkg("HMPTrees")` and `r pkg("GUniFrac")` provide tools for comparing microbial communities. `r pkg("betapart")` allows computing pair-wise dissimilarities (distance matrices) and multiple-site dissimilarities, separating the turnover and nestedness-resultant components of taxonomic (incidence and abundance based), functional and phylogenetic beta diversity. `r pkg("phyloregion")` extends  `r pkg("betapart")` to allow sparse community matrices allowing larger datasets.`r pkg("adiv")` can calculate various indices of biodiversity including species, functional and phylogenetic diversity, as well as alpha, beta, and gamma diversities. `r pkg("entropart")` can measure and partition diversity based on Tsallis entropy as well as calculate alpha, beta, and gamma diversities. `r pkg("metacoder")` is an R package for handling large taxonomic data sets, like those generated from modern high-throughput sequencing, like metabarcoding.

### Phyloclimatic modeling

`r pkg("phyloclim")` integrates several new tools in this area.

### Phylogeography and biogeography

`r pkg("phyloland")` implements a model of space colonization mapped on a phylogeny, it aims at estimating limited dispersal and competitive exclusion in a statistical phylogeographic framework. `r pkg("diversitree")` implements the GeoSSE method for diversification analyses based on two areas. `r github("GabrielNakamura/Herodotools")` can be used to perform a wide variety of biogeographical macroevolutionary analyses, including exploring spatial biodiversity patterns, conducting ancestral area reconstruction, and performing evoregion classification.

### Epidemiology

See the `r view("Epidemiology")` task view for details about packages useful for epidemiology, including phylogenetic epidemiology.

### Genetics

#### Species and population delimitation

`r pkg("adhoc")` can estimate an ad hoc distance threshold for a reference library of DNA barcodes.

#### Gene tree - species tree

`r pkg("HyPhy")` can count the duplication and loss cost to reconcile a gene tree to a species tree. It can also sample histories of gene trees from within family trees.

## Other useful packages and miscellany

### Taxonomy

`r pkg("taxize")` can interact with a suite of web APIs for taxonomic tasks, such as verifying species names, getting taxonomic hierarchies, and verifying name spelling. `r pkg("evobiR")` contains functions for making a tree at higher taxonomic levels, downloading a taxonomy tree from NCBI or ITIS, and various other miscellaneous functions (simulations of character evolution, calculating D-statistics, etc.).

### Interactions with other programs

`r pkg("geiger")` can call PATHd8 through its congruify function. `r pkg("ips")` wraps several tree inference and other programs, including MrBayes, Beast, and RAxML, allowing their easy use from within R. `r pkg("Rphylip")` wraps [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) , a broad variety of programs for tree inference under parsimony, likelihood, and distance, bootstrapping, character evolution, and more. `r pkg("BoSSA")` can use information from various tools to place a query sequence into a reference tree. `r pkg("pastis")` can use taxonomic information to make constraints for MrBayes tree searches.

### Note

 At least ten packages start as phy\* in this domain, including two pairs of similarly named packages (phytools and phylotools, phylobase and phybase). This can easily lead to confusion, and future package authors are encouraged to consider such overlaps when naming packages. For clarification, `r pkg("phytools")` provides a wide array of functions, especially for comparative methods, and is maintained by Liam Revell; `r pkg("phylotools")` has functions for building supermatrices and is maintained by Jinlong Zhang. `r pkg("phylobase")` implements S4 classes for phylogenetic trees and associated data and is maintained by Francois Michonneau; [phybase](https://code.google.com/p/phybase/) has tree utility functions and many functions for gene tree - species tree questions and is authored by Liang Liu, but no longer appears on CRAN.

## References

-   Borregaard, M.K., Rahbek, C., Fjeldsaa, J., Parra, J.L., Whittaker, R.J. and Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution 5(11): 1225-1235.
-   Butler MA, King AA 2004 Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164, 683-695.
-   Cheverud JM, Dow MM, Leutenegger W 1985 The quantitative assessment of phylogenetic constraints in comparative analyses: Sexual dimorphism in body weight among primates. Evolution 39, 1335-1351.
-   FitzJohn RG, Maddison WP, and Otto SP 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58: 595-611.
-   Garland T, Harvey PH, Ives AR 1992 Procedures for the analysis of comparative data using phylogenetically independent contrasts. Systematic Biology 41, 18-32.
-   Hansen TF 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution 51: 1341-1351.
-   Maddison WP, Midford PE, and Otto SP 2007. Estimating a binary character's effect on speciation and extinction. Systematic Biology 56: 701--710.
-   Magallon S, Sanderson, M.J. 2001. Absolute Diversification Rates in Angiosperm Clades. Evolution 55(9):1762-1780.
-   Moore, BR, Chan, KMA, Donoghue, MJ (2004) Detecting diversification rate variation in supertrees. In Bininda-Emonds ORP (ed) Phylogenetic Supertrees: Combining Information to Reveal the Tree of Life, Kluwer Academic pgs 487-533.
-   Nee S, May RM, Harvey PH 1994. The reconstructed evolutionary process. Philosophical Transactions of the Royal Society of London Series B Biological Sciences 344: 305-311.
-   Pagel M 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884
-   Pybus OG, Harvey PH 2000. Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London Series B Biological Sciences 267, 2267-2272.

### Links
-   [PHYLOCH, LAGOPUS, and PHYLOCLIM packages](http://www.christophheibl.de/Rpackages.html)
