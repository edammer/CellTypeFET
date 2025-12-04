# CellTypeFET
Calculate and Plot Significance of Gene List Overlap

Requires R packages: WGCNA, biomaRt, and RColorBrewer.

Wrapper code for the function is provided below, with examples that generate the output found in
<a href="https://github.com/edammer/CellTypeFET/blob/main/SampleOutputs.zip">SampleOutputs.zip</a>, which are
provided in this repository.

```
###################################################################################################################
# Cross-species FET modified to optionally adjust for symbol lookup inefficiency/loss (Wrapper R Script)
# by Eric Dammer & Divya Nandakumar
#=================================================================================================================#
# As published in: 
# Seyfried et al Cell Syst, 2017	https://www.sciencedirect.com/science/article/pii/S2405471216303702
# Johnson ECB et al, Nat Med, 2020	https://www.nature.com/articles/s41591-020-0815-6
# Johnson ECB et al, Nat Neurosci, 2022	https://www.nature.com/articles/s41593-021-00999-y
# ...and others
#=================================================================================================================#
# Provides and charts statistics for significance of hypergeometric overlap of gene symbol or 'UniqueID' lists.
# Can be used for cell type marker enrichment in modules, or any gene list overlap.
#
# UniqueIDs are formatted as "Symbol|ID...", where symbol is the official gene symbol (case is species specific).
#
# Sample marker list data provided are for 5 cell types of mammalian brain
# derived from thresholded filtering of supplemental data in acutely
# isolated purified cells from mouse brain measured as mRNA, in:
#  Ye Zhang et al, J Neurosci, 2014
#  https://www.jneurosci.org/content/34/36/11929.short
#  (Barres RNA app now at https://www.brainrnaseq.org/ )
#
# or measured as protein, in:
#  Kirti Sharma et al, Nat Neurosci, 2015 https://www.nature.com/articles/nn.4160
#
# Sample gene clusters are provided for single cell cluster markers in mRNA
# from fly brain, as published in:
#  Kristofer Davie et al, Cell, 2018
#  https://www.sciencedirect.com/science/article/pii/S0092867418307207
#
# Sample WGCNA modules representing proteome organization of the motor cortex
# of the brains of ALS and healthy individuals is provided for testing use only.
# Please reach out if you are interested in repurposing this unpublished data.
#=================================================================================================================#
# edammer@emory.edu Eric Dammer, Bioinformatic Scientist, NT Seyfried Systems Biology Group
# Emory University School of Medicine (2023)
#=================================================================================================================#


############ PREPARATION OF IN-MEMORY VARIABLES (optional) ##########################
rootdir="./"  				# Contains inputs and outputs
setwd(rootdir)

## Load into memory 3 variables for your network
load("motorCtx-ALS_data_forORA.Rdata")  # sample Rdata file provided to demonstrate the script using motor cortex brain total homogenate in an unpublished ALS cohort
					                    # DO NOT RELEASE OR USE WITHOUT PERMISSION

## Standardize the variable names of for modulesInMemory option below
net<-netALS  		    # list output of WGCNA::blockwiseModules() function building your coexpression network. The one necessary item in the list is the colors variable, i.e. net$colors or net[["colors"]]
cleanDat<-cleanDatALS	# data.frame with rows representing gene products measured and all making it into the WGCNA network as goodGenes.
			            # Rownames are expected to have the official gene symbol followed by a pipe and any other information after, i.e. "NEFL|..."
numericMeta<-numericMetaALS   #sample traits/metadata


#################### CONFIGURATION PARAMETERS FULL LIST #############################
##             WITH SAMPLE VALUES GIVEN FOR SAMPLE DATA PROVIDED                   ##

heatmapScale="minusLogFDR"                        				# Accepted options are "p.unadj" or "minusLogFDR"
heatmapTitle="My Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists"	# What are your categories (or WGCNA) list of lists based on?
                                                                                        # And What gene lists are your reference lists?
paletteColors="YlGnBu"                                          # See valid palettes using RColorBrewer::display.brewer.all()
                                                                # Can be a vector if there are more than 1 refDataFiles (heatmaps to generate)

FileBaseName="MyNetworkModules_FET_to_5brainCellTypeMarkerLists"
refDataDescription="5brainCellTypes"				# One Description of reference Data list(s) specified in PDF file name below
# File Names of Reference List(s): You will get one output PDF page per file
refDataFiles <- c(      "MyGene-Human-SharmaZhangUnion.csv",	# HUMAN gene symbols have been pre-converted from the original below MOUSE lists.
			            "MyGene-Mouse-SharmaZhangUnion.csv")	# Originally mined brain cell type mRNA and protein lists were based on experiments in MOUSE.
speciesCode=c("hsapiens","mmusculus") 				# species code(s) for biomaRt (one for each refDataFile)#one for each .csv in refDataFiles


# Use Modules in memory OR a .csv file with your input gene lists.
modulesInMemory=TRUE                              		    # Load modules as categories? (If TRUE, categoriesFile not used, but you need cleanDat, net[["colors"]] and numericMeta variables)
categoriesFile="Fly_Seurat_87cluster_BrainCellTypes.csv"	# File Name of Categories (Lists of Fly genes), only loaded if modulesInMemory=FALSE
								                            # NOTE this file format has a column for official gene symbols of each module or cluster, with the cluster name/ID as column names in row 1
categorySpeciesCode="dmelanogaster"				            # What species are the gene sybmols in categoriesFile?
bkgrFileForCategories="Dmel_geneBackground.csv"             # 12/2025: Augment background of categoriesFile symbols with a single-column no-header list (in a no comma .CSV file)

# Other Options
allowDuplicates=TRUE				# Allow duplicate symbols across different row lists (refDataFiles) for overlap?
                                    # if FALSE, row lists' duplicte entries removed entirely (since removing only extra instance(s) would have to choose to keep in one row)
						            # (should be true if you have general cell type lists and e.g. disease-associated phenotype cell type lists)
strictSymmetry=FALSE                # remove duplicated symbols across columns (modules or categoriesFile, background lists)
                                    # if TRUE, and all row list symbols are present in column lists, background, then swapping row and column input will give similar or identical stats
resortListsDecreasingSize=FALSE		# resort categories/modules and reference data lists? (decreasing size order)
heatmapScale="FDR"                  # "FDR" (or still accepted previous default "minusLogFDR") will FDR-correct all p values used for plotting. Otherwise "p" or "p.unlog". This is independent of the legendScale setting.
legendScale="minusLog"              # Heatmap legend will be plotted with increasing significance and color intensity for higher -log10(p or FDR). "unlog" option has increasing significance for lower p or FDR.
barOption=FALSE					    # Draw bar charts for each list overlap instead of a heatmap.
asterisksOnly=FALSE                 # 12/2025: new (TRUE, default); if FALSE, print p (or FDR) values within heatmap cells if below maxPcolor and expand width of columns to accomodate this text
maxPcolor=0.25                      # 12/2025: p (or FDR) values below this value will get some heatmap color in the scale; above the value, cells will be WHITE. Plots with no values below this value only show color in the legend.
adjustFETforLookupEfficiency=FALSE	# adjust p FET input for cross-species lookup inefficiency/loss of list member counts?
verticalCompression=1				# DEPRECATED for values >1. Plot(s) may be squeezed into 1 row out of this many in each PDF page, compressing the heatmap tracks vertically (or the bar chart heights) for each reference list)
                                    # (PDF dimensions and margins for plots are now dynamically calculated to fit all plot elements and text optimally, one per page.)
reproduceHistoricCalc=FALSE			# should be FALSE unless trying to reproduce exact calculations of prior publications listed.
#####################################################################################


## Generate Sample Outputs

# Load Seyfried/Emory pipeline FET as function geneListFET() having all the parameters described above, many with defaults used.
source("./geneListFET.R")


## output enrichment significance of provided test network data in memory as -log10(FDR) heatmap; enrichments checked are in provided 5 brain cell type marker gene lists
geneListFET(FileBaseName="1.ALSnetModules_FET_to_5brainCellTypes",
            heatmapTitle="ALS Motor Cortex Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists",
            maxPcolor=maxPcolor, asterisksOnly=asterisksOnly, paletteColors="plasma",  # palettes from viridisLite package valid, as well as RcolorBewer
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode=c("hsapiens","mmusculus"),refDataDescription="5brainCellTypes")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?


## Same as above, but output a bar chart instead of a heatmap.
geneListFET(FileBaseName="2.ALSnetModules_FET_to_5brainCellTypes.barChart", barOption=TRUE,
            heatmapTitle="ALS Motor Cortex Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists",
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  # use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode=c("hsapiens","mmusculus"),refDataDescription="5brainCellTypes")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?
#NOTE: Sample output #2 has 5 bar charts for the ALS motor cortex network's modules' gene enrichment into each of the human of 5 cell type marker lists.
#      The next 5 pages have the same network's modules' gene enrichment in the mouse symbol lists of the same 5 cell type markers. Total PDF page count is 10.


## output enrichment significance of provided "categories" (or fly brain gene clusters) input from a .csv file instead of using coexpression modules in memory; -log10(FDR) heatmap; enrichments are checked against provided 5 mammalian brain cell type marker gene lists
geneListFET(FileBaseName="3.FlyBrainSCgene87SeuratClusters_FET_to_5mammalianBrainCellTypes-symmetry",
            heatmapTitle="Fly Brain 87 Cluster Overlaps with 5 Mammalian Brain Cell Type Marker Reference Lists",
            maxPcolor=maxPcolor, asterisksOnly=asterisksOnly, paletteColors="plasma",  # palettes from viridisLite package valid, as well as RcolorBewer
            strictSymmetry=TRUE,  # When TRUE, categories (lists for columns) duplicate genes are only considered once (dereplicated to one overlap per gene in the backaground)
            modulesInMemory=FALSE,categoriesFile="Fly_Seurat_87cluster_BrainCellTypes.csv",categorySpeciesCode="dmelanogaster",  #use clusters of gene symbols in a provided file, categoriesFile; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            bkgrFileForCategories="Dmel_geneBackground.csv",  # can have a large effect on significance range, by augmenting the background or universe for overlap (in this case, to protein-coding fly symbols)
            refDataFiles=refDataFiles,speciesCode=c("hsapiens","mmusculus"),refDataDescription="5brainCellTypes")  #file(s) with columns of reference gene lists to check for overlap in; what are the species code(s) for symbols in each file?


## Duplicated genes in >1 reference list are not explicitly checked for, nor removed, currently. Instances occur due to different gene splicing isoforms specific to different cell types, e.g.
#  We provide our human reference list in 3 versions: (1) the reference with maximal duplication n=132 duplicated symbols provided on GitHub originally;
#                                                     (2) the published reference 5 cell type human list re-converted from mouse symbols using mygene for publication in 2022;
#                                                     (3) the reference 5 cell types from #2 with the 19 genes duplicated in those removed.
# These lists are compared (with heatmap output):
geneListFET(FileBaseName="1.ALSnet_dupCtl",
            heatmapTitle="ALS Motor Cortex Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists",
            maxPcolor=maxPcolor, asterisksOnly=asterisksOnly, paletteColors="plasma",  # palettes from viridisLite package valid, as well as RcolorBewer
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",
            refDataFiles=c("MyGene-Human-SharmaZhangUnion.csv",
                           "MyGene-Human-SharmaZhangUnion_Published(JohnsonECB_NatNeurosci_2022).csv",
                           "MyGene-Human-SharmaZhangUnion_Published(JohnsonECB_NatNeurosci_2022)_uniqueONLY.csv"),
            speciesCode=c("hsapiens","hsapiens","hsapiens"),refDataDescription="5brainCellTypes")
```
## Note on duplicate gene symbols in _different_ reference cell type marker lists

Exact reference lists from Johnson ECB, et al (2022) are now provided as converted from mouse using symbol interspecies lookup using the mygene package, and a version of the lists not published with 19 duplicated genes removed '(uniqueONLY)' is also added to this repository. The above check shows similar significance for enrichments across all 3 of the lists, due to the generally unchanged overlaps, overall.

## Notes on handling of duplicate gene symbols in general

From experience, it is important to be aware of and handle gene symbol redundancy within and across lists (both row and column) to get expected, consistent, and statistically sound results.

Reference lists (rows, or barplot inputs per barplot page) input as .csv files do NOT contribute to background. Genes in these lists not overlapping with category lists (symbols in lists representing columns or bars) do not impact the contingency table counts for fisher.test. If you want to see the same p values or FDR values with reference (row) and category (column) lists swapped, be sure that all reference list genes are in the category lists (e.g. grey gene list not assigned to a module), and handle duplicates strictly. When modulesInMemory=FALSE and a categoriesFile is specified for column-specific lists, an additional file for augmentation or completion of the gene 'universe' or background can be supplied via the **bkgrFileForCategories** argument. To remove all redundancy (duplicate gene symbols) across all categories and any category background, use the argument **strictSymmetry=TRUE**. This is not recommended for coexpression module-derived categories/columns, as multiple isoforms for the same gene in different modules will be undercounted, but is recommended if your column-specific lists are highly redundant--e.g., representing an evolving time-window series of gene markers. To remove duplicates within (but not across) reference lists--which are considered as independent marker lists--use the argument **allowDuplicates=FALSE**.

## Package Dependencies

WGCNA, biomaRt, RColorBrewer, viridisLite (either of the last 2 is required, depending on which palette is named for argument **paletteColors**; default is "YlGnBu" palette of RColorBrewer).
