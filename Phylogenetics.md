---
name: Phylogenetics
topic: Phylogenetics packages in R
maintainer: William Gearty, Brian O'Meara, Jacob Berv, Gustavo Ballen Chaparro, Diniz Ferreira, Hilmar Lapp, Lars Schmitz, Martin R. Smith, Nathan S. Upham
email: willgearty@gmail.com
version: 2022-09-14
source: https://github.com/cran-task-views/Phylogenetics/
---

# Overview

The history of life unfolds within a phylogenetic context, and phylogenetic trees (often shortened to "trees") are developed to represent this evolutionary history. Comparative phylogenetic methods are statistical approaches for analyzing historical patterns along such phylogenetic trees. This task view describes the plethora of R packages that 1) are useful for the handling and manipulation of phylogenetic trees and/or 2) implement a variety of different comparative phylogenetic methods. This is an active research area and much of the information is subject to change. One thing to note is that many important packages are not on CRAN: either they were formerly on CRAN and were later archived (for example, if they failed to incorporate necessary changes as R is updated) or they are developed elsewhere and have not been put on CRAN yet. Such packages may be found on GitHub, R-Forge, [Bioconductor](https://bioconductor.org/packages/release/BiocViews.html#___Phylogenetics), or authors' websites. Another important note is that at least ten packages start as phy\* in this domain, including two pairs of similarly named packages (phytools and phylotools, phylobase and phybase). Users are encouraged to read and distinguish carefully between the names of the packages that are listed below.

If you have any questions, feel free to reach out to the task view maintainers or the maintainers of specific packages. Questions may also be directed to the [R-SIG-Phylo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Phylo/) mailing-list after subscription.

# Core packages

- `r pkg("ape", priority = "core")` implements the S3 phylo class which is commonly used to store phylogenetic trees in R. It has numerous methods for reading, writing, and plotting trees; manipulating and building trees in R; and performing phylogenetic comparative analyses. These methods are used as a framework for many other phylogenetics R packages.
- `r pkg("phylobase", priority = "core")` implements the S4 phylo4 class which combines phylogenetic trees and comparative data. While not used as commonly as the S3 phylo class, this new class is gaining traction among newer packages that implement phylogenetic comparative methods (e.g., `r pkg("adephylo")` and `r pkg("phylosignal")`).
- `r pkg("geiger", priority = "core")` implements a large suite of model fitting approaches for analyses of trait evolution and diversification. It also has several auxiliary functions that are often used by other packages.
- `r pkg("phytools", priority = "core")` has a constantly increasing range of functions for performing phylogenetic comparative methods and visualizing, manipulating, reading or writing, and even inferring phylogenetic trees.

# Other packages

Non-core packages are grouped into the following categories:

1.  **Working with trees in R:** packages dedicated to the handling, manipulation, and visualization of phylogenetic data
2.  **Building trees in R:** packages for phylogenetic inference and tree simulation
3.  **Comparative phylogenetic methods:** packages for performing various comparative phylogenetic methods, including those dealing with trait evolution and diversification
4.  **Phylogenetics in other fields:** packages designed to perform field-specific phylogenetic analyses, including paleontology, community ecology, biogeography, and genetics
4.  **Other useful packages and miscellany:** packages that are useful for performing phylogenetic analyses, such as taxonomic matching

## Working with trees in R

### Getting trees into R

- `r pkg("ape")` can read trees from external files in newick format (sometimes popularly known as phylip format) or NEXUS format. It can also read trees input by hand as a newick string (e.g., "(human,(chimp,bonobo));").
- `r pkg("phylobase")` and its lighter weight sibling `r pkg("rncl")` can use the [Nexus Class Library](http://ncl.sourceforge.net/) to read NEXUS, Newick, and other tree formats.
- `r pkg("treebase")` can search for and load trees from the online tree repository TreeBASE.
- `r pkg("rdryad")` can pull data from the online data repository Dryad.
- `r pkg("RNeXML")` can read, write, and process metadata for the [NeXML](http://www.nexml.org) format.
- `r pkg("TreeTools")` can read trees from external files in TNT format and NEXUS format, including extensions to the Nexus format not supported by ape, and metadata from MorphoBank.
- PHYLOCH can load trees from BEAST, MrBayes, and other phylogenetics programs (PHYLOCH is only available from the author's [website](http://www.christophheibl.de/Rpackages.html) ).
- `r pkg("phyext2")` can read and write various tree formats, including simmap formats.
- `r pkg("rotl")` can pull in a synthetic tree and individual study trees from the Open Tree of Life project.
- The `r bioc("treeio")` package can read trees in Newick, Nexus, New Hampshire eXtended format (NHX), jplace and Phylip formats and data output from BEAST, EPA, HyPhy, MrBayes, PAML, PHYLDOG, pplacer, r8s, RAxML and RevBayes.
- `r pkg("phylogram")` can convert Newick files into dendrogram objects (see `r pkg("dendrogram")` for the manipulation of such objects).
- `r pkg("brranching")` can fetch phylogenies from online repositories, including [phylomatic](http://phylodiversity.net/phylomatic/).

### Tree manipulation

- `r pkg("ape")` has functions for rooting trees, dropping tips, randomly resolving polytomies, creating branch lengths, getting information about tree size or other properties, pulling in data from GenBank, and many more.
- `r pkg("phylobase")` has functions for traversing a tree (e.g., getting all descendants from a particular node specified by just two of its descendants).
- `r pkg("geiger")` can prune trees and data to an overlapping set of taxa. It can be also used to perform branch length scaling using ACDC; Pagel's (1999) lambda, delta and kappa parameters; and the Ornstein-Uhlenbeck alpha parameter (for ultrametric trees only). It can also be used to prune extinct taxa.
- `r pkg("TreeTools")` has functions to quantify and manipulate tree shape and balance, including the application of constraints; and to measure the phylogenetic information content of trees. Rogue identifies wildcard taxa, generating more informative summary trees.
- `r pkg("tidytree")` can convert a tree object to a tidy data frame and has other tidy approaches to manipulate tree data.
- `r pkg("evobiR")` can do fuzzy matching of names (to allow some differences).
- `r pkg("SigTree")` finds branches that are responsive to some treatment, while allowing correction for multiple comparisons.
- `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more.
- `r pkg("apex")` can handle multiple gene DNA alignments making their use and analysis for tree inference easier in `r pkg("ape")` and `r pkg("phangorn")`.
- `r pkg("aphid")` can weight sequences based on a phylogeny and can use hidden Markov models (HMMs) for a variety of purposes including multiple sequence alignment.
- `r pkg("phytools")` also allows branch length scaling, as well as several tree transformations (adding tips, finding subtrees).
- `r pkg("phangorn")` and `r pkg("TreeSearch")` can perform tree rearrangements (NNI, SPR, and TBR).
- `r pkg("paleotree")` has functions for manipulating trees based on sampling issues that arise with fossil taxa as well as more universal transformations.
- `r pkg("dendextend")` can manipulate dendrograms, including subdividing trees, adding leaves, and more.

### Tree visualization

- `r pkg("ape")`, `r pkg("adephylo")`, `r pkg("phylobase")`, `r pkg("phytools")`, `r pkg("ouch")`, and `r pkg("dendextend")` have functions for plotting trees; several of these have options for branch or taxon coloring based on some criterion (ancestral state, tree structure, etc.). Trees can be examined (zoomed) and viewed as correlograms using `r pkg("ape")`.
- `r rforge("paleoPhylo")` and `r pkg("paleotree")` are specialized for drawing paleobiological phylogenies.
- `r pkg("phytools")` can project a tree into a morphospace.
- The popular R visualization package `r pkg("ggplot2")` can be extended by `r bioc("ggtree")` and `r bioc("ggtreeExtra")` to visualize phylogenies, and a geological timescale can be added using `r pkg("deeptime")`.
- `r pkg("strap")` can be used to add a geological timescale to a phylogeny, along with stratigraphic ranges.
- `r pkg("idendr0")` can be used to interactively explore trees (as dendrograms).
- `r pkg("phylocanvas")` is a widget for "htmlwidgets" that enables embedding of phylogenetic trees using the phylocanvas javascript library.
- `r pkg("ggmuller")` allows plotting a phylogeny along with frequency dynamics.
- `r pkg("RPANDA")` can be used to plot the spectral density and eigenvalues of a phylogeny.

### Tree comparison

- `r pkg("distory")`, `r pkg("Quartet")`, `r pkg("TBRDist")` and `r pkg("TreeDist")` can be used to evaluate tree-tree distances.
- `r pkg("ape")` can compute tree-tree distances and also create a plot showing two trees with links between associated tips.
- `r pkg("kdetrees")` implements a non-parametric method for identifying potential outlying observations in a collection of phylogenetic trees, which could represent inference problems or processes such as horizontal gene transfer.
- `r pkg("dendextend")` can evaluate multiple measures comparing dendrograms.

## Tree building in R

### Phylogenetic inference

- `r pkg("ape")` can be used to perform neighbour joining, bio-nj and fast ME methods of phylogenetic reconstruction.
- `r pkg("phangorn")` can estimate trees using distance (e.g. UPGMA), parsimony, and likelihood.
- `r pkg("TreeSearch")` can identify most-parsimonious trees under parsimony, including with inapplicable data, and includes a graphical user interface for detailed analysis of results.
- `r pkg("phyclust")` can cluster sequences.
- `r pkg("phytools")` can build trees using MRP supertree estimation and least squares.
- `r pkg("phylotools")` can build supermatrices for analyses in other software.
- `r pkg("pastis")` can use taxonomic information to make constraints for Bayesian tree searches.
- `r pkg("EvoPhylo")` can be used to perform automated morphological character partitioning for bayesian phylogenetic analyses that are performed with [MrBayes](http://nbisweden.github.io/MrBayes/) and [BEAST2](https://www.beast2.org/). It can also be used to analyze the macroevolutionary parameter outputs from such analyses.
- `r bioc("fastreeR")` can be used to calculate distances, build phylogenetic trees, or perform hierarchical clustering between the samples of a VCF or FASTA file.

### Divergence times

- `r pkg("ape")` implements non-parametric rate smoothing (NPRS) and penalized likelihood.
- `r pkg("geiger")` can do congruification to stretch a source tree to match a specified standard tree.
- `r pkg("treedater")` implements various clock models, ways to assess confidence, and detecting outliers.
- `r pkg("phangorn")` can infer ultrametric and tipdated phylogenies with a strict clock model direct from sequences.  

### Tree simulations

- `r pkg("TreeSim")` can be used to simulate trees using constant-rate birth-death with various constraints.
- `r pkg("geiger")` can be used to simulate trees under a birth-death process.
- `r pkg("ape")` can be used to generate random trees by random splitting of edges (for non-parametric trees) or random clustering of tips (for coalescent trees).
- `r pkg("paleotree")` can simulate fossil deposition, sampling, and the tree arising from this as well as trees conditioned on observed fossil taxa.
- `r pkg("FossilSim")` can be used to simulate fossil data on existing phylogenetic trees under mechanistic models of preservation and sampling.
- `r pkg("TESS")` can simulate trees with time-dependent speciation and/or extinction rates, including mass extinctions.
- `r pkg("paleobuddy")` presents a flexible interface to simulate a wide array of user-defined diversification dynamics, including environmental-dependence.

## Comparative phylogenetic methods

### Ancestral state reconstruction

- `r pkg("ape")` can reconstructe continuous characters using maximum likelihood, generalised least squares, or independent contrasts.
- `r pkg("ouch")` can be used to reconstruct root ancestral character states under Brownian motion or Ornstein-Uhlenbeck models, though ancestral states at the internal nodes are not.
- `r pkg("ape")` can reconstruct discrete characters using a variety of Markovian models that parameterize the transition rates among states.
- `r pkg("markophylo")` can fit a broad set of discrete character types with models that can incorporate constrained substitution rates, rate partitioning across sites, branch-specific rates, sampling bias, and non-stationary root probabilities.
- `r pkg("phytools")` can do stochastic character mapping of traits on trees.
- `r pkg("Rphylopars")` can perform ancestral state reconstruction for datasets with multiple observations per species and/or missing data.
- `r pkg("TreeSearch")` can perform mapping of characters under parsimony, with an allowance for inapplicable data.

### Trait evolution

- `r pkg("ape")`, `r pkg("picante")`, or `r pkg("caper")` can be used to calculate independent contrasts for continuous characters. `r pkg("caper")` also implements the brunch and crunch algorithms.
- `r pkg("geiger")` can be used to perform analyses of discrete trait evolution, including models of unequal rates or rates changing at a given instant of time, as well as Pagel's transformations.
- `r pkg("geiger")`, `r pkg("ape")`, `r pkg("paleotree")`, and `r pkg("motmot")` can be used to fit Brownian motion models.
- `r github("cran/RBrownie")` can fit multiple-rate Brownian motion models.
- `r pkg("geiger")` and `r pkg("OUwie")` can be used to investigate deviations from Brownian motion.
- `r pkg("geiger")` can be used to fit other continuous models, including Pagel's transforms and models with trends.
- `r pkg("mvMORPH")` can fit Brownian motion, early burst, ACDC, OU, and shift models to univariate or multivariate data.
- `r pkg("geiger")`, `r pkg("ape")`, `r pkg("motmot")`, `r pkg("ouch")`, `r pkg("surface")`, and `r pkg("OUwie")` can be used to fit Ornstein-Uhlenbeck (OU) models. `r pkg("ouch")` can implement models with multiple means, `r pkg("surface")` can implement models with multiple means using stepwise AIC, and `r pkg("OUwie")` can implement models with multiple means, rates, and attraction values. Also see `r github("mongiardino/extendedSurface")` which combines the functionality of `r pkg("surface")` and `r pkg("OUwie")`.
- `r pkg("motmot")` can be used to fit continuous models that change rate or mode at specific time(s).
- `r pkg("Rphylopars")` can be used to fit continuous models such as those described above to datasets with multiple observations per species and/or missing data.
- `r pkg("geiger")` implements ANOVA's and MANOVA's in a phylogenetic context.
- `r pkg("ape")`, `r pkg("PHYLOGR")`, `r pkg("caper")`, and `r pkg("motmot")` implement traditional GLS methods (sensu Grafen or Martins).
- `r pkg("ape")` can be used to calculate phylogenetic autoregression (sensu Cheverud et al).
- `r pkg("ape")` and `r pkg("adephylo")` can be used to calculate phylogenetic autocorrelation (Moran's I).
- `r pkg("MCMCglmm")` can be used to assess correlation between traits using a GLMM.
- `r pkg("phylolm")` can fit phylogenetic linear regression and phylogenetic logistic regression models using a fast algorithm, making it suitable for large trees.
- `r pkg("brms")` can examine correlations between continuous and discrete traits, and can incorporate multiple measurements per species.
- `r pkg("phytools")` can also investigate rates of trait evolution and do stochastic character mapping.
- `r pkg("metafor")` can perform meta-analyses accounting for phylogenetic structure.
- `r pkg("pmc")` evaluates the model adequacy of several trait models (from `r pkg("geiger")` and `r pkg("ouch")`) using Monte Carlo approaches.
- `r pkg("phyreg")` implements the Grafen (1989) phyglogenetic regression.
- `r pkg("geomorph")` can do geometric morphometric analysis in a phylogenetic context.
- `r pkg("dispRity")` can be used to calculate disparity through time and perform other disparity-related analyses.
- `r pkg("MPSEM")` can predict features of one species based on information from related species using phylogenetic eigenvector maps.
- `r pkg("convevol")` and `r pkg("windex")` can both test for convergent evolution on a phylogeny.
- `r pkg("Claddis")` can be used to measure morphological diversity from discrete character data and evolutionary tempo on phylogenetic trees.

### Trait simulations

- `r pkg("ouch")`, `r pkg("geiger")`, `r pkg("ape")`, `r pkg("picante")`, `r pkg("OUwie")`, and `r pkg("caper")` can be used to simulate continuous traits using Brownian motion.
- `r pkg("ouch")` and `r pkg("OUwie")` can be used to simulate continuous traits using the Hansen model (a form of the OU).
- `r pkg("geiger")` can be used to simulate continuous traits using a speciational model and discrete traits can be simulated using a continuous time Markov model (including models where rates change through time).
- `r pkg("phangorn")` can simulate DNA or amino acids.
- `r pkg("phytools")` can simulate discrete characters using stochastic character mapping.
- `r pkg("phylolm")` can simulate continuous or binary traits along a tree.
- `r pkg("Rphylopars")` can simulate data with missing observations.

### Diversification analysis

- `r pkg("ape")` can be used to create lineage through time plots.
- `r pkg("ape")` can fit a simple birth-death model for when you have extant species only (sensu Nee et al. 1994), survival models, and goodness-of-fit tests (as applied to testing of models of diversification).
- `r pkg("TESS")` can calculate the likelihood of a tree under a model with time-dependent diversification, including mass extinctions.
- `r pkg("geiger")` can calculate net rates of diversification (sensu Magellon and Sanderson).
- `r pkg("diversitree")` implements the BiSSE method (Maddison et al. 1997) and later improvements (FitzJohn et al. 2009).
- `r pkg("hisse")` implements various hidden state diversification models, including HiSSE (Beaulieu and O'Meara 2016), GeoHiSSE (Caetano et al. 2018), MuHiSSE (Nakov et al. 2019), and MiSSE (trait-independent).
- `r pkg("TreePar")` estimates speciation and extinction rates with models where rates can change as a function of time (i.e., at mass extinction events) or as a function of the number of species.
- `r pkg("caper")` can do the macrocaic test to evaluate the effect of a a trait on diversity.
- `r pkg("apTreeshape")` also has tests for differential diversification (see [description](https://doi.org/10.1093/bioinformatics/bti798) ).
- `r pkg("iteRates")` can identify and visualize areas on a tree undergoing differential diversification.
- `r pkg("DDD")` implements maximum likelihood methods based on the diversity-dependent birth-death process to test whether speciation or extinction are diversity-dependent, as well as identifies key innovations and simulate a density-dependent process.
- `r pkg("PBD")` can calculate the likelihood of a tree under a protracted speciation model.
- `r pkg("phyloTop")` has functions for investigating tree shape, with special functions and datasets relating to trees of infectious diseases.
- `r pkg("RPANDA")` can be used to fit various diversification models to phylogenies, including time-dependent and environmental-dependent models.

## Phylogenetics in specific fields

### Morphometrics

- `r pkg("geomorph")` and `r pkg("RRPP") may be used to evaluate evolutionary trends in multivariate phenotypes. Available methods include phylogenetic linear models (phylogenetic anova/regression, etc.), phylogenetic partial least squares, comparing rates of phenotypic evolution, phylogenetic integration, and phylogenetic modularity. Additionally, ordination approaches include both phylogenetic PCA, and phylogenetically-aligned components analysis (PACA).

### Time series and paleontology

- `r pkg("paleoTS")` can be used to analyze paleontological time series data using a likelihood-based framework for fitting and comparing models (using a model testing approach) of phyletic evolution (based on the random walk or stasis model).
- `r pkg("strap")` can do stratigraphic analysis of phylogenetic trees.
- `r github("rachelwarnock/fbdR")` can be used to estimate diversification rates from phylogenetic trees and fossil occurrence data.
- R offers a wealth of other options for general-purpose time series modeling, many of which are listed in the `r view("TimeSeries")` task view.

### Community and microbial ecology

- `r pkg("picante")`, `r pkg("vegan")`, `r pkg("SYNCSA")`, `r pkg("phylotools")`, `r pkg("PCPS")`, `r pkg("caper")`, `r pkg("DAMOCLES")`, and `r pkg("phyloregion")` integrate several tools for using phylogenetics with community ecology.
- `r pkg("HMPTrees")` and `r pkg("GUniFrac")` provide tools for comparing microbial communities.
- `r pkg("betapart")` allows computing pair-wise dissimilarities (distance matrices) and multiple-site dissimilarities, separating the turnover and nestedness-resultant components of taxonomic (incidence and abundance based), functional and phylogenetic beta diversity.
- `r pkg("phyloregion")` extends  `r pkg("betapart")` to allow sparse community matrices allowing larger datasets.
- `r pkg("adiv")` can calculate various indices of biodiversity including species, functional and phylogenetic diversity, as well as alpha, beta, and gamma diversities.
- `r pkg("entropart")` can measure and partition diversity based on Tsallis entropy as well as calculate alpha, beta, and gamma diversities.
- `r pkg("metacoder")` provides functions for handling large taxonomic data sets, like those generated from modern high-throughput sequencing, like metabarcoding.
- `r bioc("phyloseq")` provides a set of classes and tools to facilitate the import, storage, analysis, and graphical display of microbiome census data.

### Phyloclimatic modeling

- `r pkg("phyloclim")` integrates several new tools in this area.

### Phylogeography and biogeography

- `r pkg("phyloland")` implements a model of space colonization mapped on a phylogeny, it aims at estimating limited dispersal and competitive exclusion in a statistical phylogeographic framework.
- `r pkg("diversitree")` implements the GeoSSE method for diversification analyses based on two areas.
- `r github("GabrielNakamura/Herodotools")` can be used to perform a wide variety of biogeographical macroevolutionary analyses, including exploring spatial biodiversity patterns, conducting ancestral area reconstruction, and performing evoregion classification.
- `r pkg("epm")` (EcoPhyloMapper) can be used to calculate various morphological and phylogenetic community metrics across geography.

### Epidemiology

See the `r view("Epidemiology")` task view for details about packages useful for epidemiology, including phylogenetic epidemiology.

### Gene tree - species tree

- `r pkg("HyPhy")` can count the duplication and loss cost to reconcile a gene tree to a species tree. It can also sample histories of gene trees from within family trees.

## Other useful packages and miscellany

### Taxonomy

- `r pkg("taxize")` can interact with a suite of web APIs for taxonomic tasks, such as verifying species names, getting taxonomic hierarchies, and verifying name spelling.
- `r pkg("evobiR")` contains functions for making a tree at higher taxonomic levels, downloading a taxonomy tree from NCBI or ITIS, and various other miscellaneous functions (simulations of character evolution, calculating D-statistics, etc.).

### Interactions with other programs

- `r pkg("ape")` can call [PhyML](http://www.atgc-montpellier.fr/phyml/), [Clustal](http://www.clustal.org/), [T-Coffee](https://tcoffee.crg.eu/), and [Muscle](https://www.drive5.com/muscle/) through various functions.
- `r pkg("geiger")` can call PATHd8 through its congruify function.
- `r pkg("ips")` wraps several tree inference and other programs, including MrBayes, Beast, and RAxML, allowing their easy use from within R.
- `r pkg("Rphylip")` wraps [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) , a broad variety of programs for tree inference under parsimony, likelihood, and distance, bootstrapping, character evolution, and more.
- `r pkg("BoSSA")` can use information from various tools to place a query sequence into a reference tree.
- `r pkg("pastis")` can use taxonomic information to make constraints for MrBayes tree searches.
- `r pkg("BAMMtools")` is an interface to the BAMM program to allow visualization of rate shifts, comparison of diversification models, and other functions.
- `r pkg("Revticulate")` can be used to interact with [RevBayes](https://revbayes.github.io/) from within R, while `r pkg("RevGadget")` can be used to process the output generated by RevBayes.

## References

-   Beaulieu, J.M. and O’Meara, B.C., 2016. Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. Systematic Biology, 65(4): 583-601.
-   Borregaard, M.K., Rahbek, C., Fjeldsaa, J., Parra, J.L., Whittaker, R.J. and Graham, C.H. 2014. Node-based analysis of species distributions. Methods in Ecology and Evolution 5(11): 1225-1235.
-   Butler MA, King AA 2004 Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164, 683-695.
-   Caetano, D.S., B.C. O’Meara, and J.M. Beaulieu. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographic models. Evolution, 72:2308-2324.
-   Cheverud JM, Dow MM, Leutenegger W 1985 The quantitative assessment of phylogenetic constraints in comparative analyses: Sexual dimorphism in body weight among primates. Evolution 39, 1335-1351.
-   FitzJohn RG, Maddison WP, and Otto SP 2009. Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. Systematic Biology 58: 595-611.
-   Garland T, Harvey PH, Ives AR 1992 Procedures for the analysis of comparative data using phylogenetically independent contrasts. Systematic Biology 41, 18-32.
-   Hansen TF 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution 51: 1341-1351.
-   Maddison WP, Midford PE, and Otto SP 2007. Estimating a binary character's effect on speciation and extinction. Systematic Biology 56: 701--710.
-   Magallon S, Sanderson, M.J. 2001. Absolute Diversification Rates in Angiosperm Clades. Evolution 55(9):1762-1780.
-   Moore, BR, Chan, KMA, Donoghue, MJ (2004) Detecting diversification rate variation in supertrees. In Bininda-Emonds ORP (ed) Phylogenetic Supertrees: Combining Information to Reveal the Tree of Life, Kluwer Academic pgs 487-533.
-   Nakov, T., Beaulieu, J.M., and Alverson, A.J. 2019. Diatoms diversify and turn over faster in freshwater than marine environments. Evolution, 73: 2497-2511.
-   Nee S, May RM, Harvey PH 1994. The reconstructed evolutionary process. Philosophical Transactions of the Royal Society of London Series B Biological Sciences 344: 305-311.
-   Pagel M 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884
-   Pybus OG, Harvey PH 2000. Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London Series B Biological Sciences 267, 2267-2272.

### Links
-   [PHYLOCH, LAGOPUS, and PHYLOCLIM packages](http://www.christophheibl.de/Rpackages.html)
