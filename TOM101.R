## lOAD REQUIRED PACKAGES
#install.packages(c())
library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library(tidyverse)
library('flashClust')
library('igraph')

# Make sure results are reproducible
set.seed(1)

# Read sample meta-data file
samples <- read.csv('/home/user/Desktop/SraRunTable.csv', header=TRUE)
samples<-samples[,c(1,25)]

#see results after retaining most important columns (sample name and serologic_response_status)
kable(samples)

names(samples)[1]<-"sample_id"
names(samples)[2]<-"condition"
kable(samples)

samples$condition<-factor(samples$condition)
levels(samples$condition)[levels(samples$condition)=="unvaccinated"] <- "unvac"
levels(samples$condition)[levels(samples$condition)=="vaccinated, not protected"] <- "vac_nt_pr"
levels(samples$condition)[levels(samples$condition)=="vaccinated, protected"] <- "vac_pr"

samples$condition

# Read in count data
raw_counts <- read.table('/home/user/Desktop/NA.txt', row.names=1, header=TRUE)
head(raw_counts)
dim(raw_counts)

names(raw_counts)
#raw_counts<-raw_counts[-c(1:5),]
#wkdir<"/home/user/Desktop/results/"
#setwd(wkdir)

## preprocessing bit. Remove low count genes
# should atleast ahve a count of 1 read per sample (atleast 45 reads per gene)

low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)

sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask),
        sum(!low_count_mask))

# Log2 transformation
log_counts <- log2(raw_counts + 1)

# Remove non differentially-expressed genes
# first, let's remove any genes with _zero_ variance since these are not
# going to help us, and may cause problems with some of the models
# remove genes that have zero variance as they will have unifoem expression over time, this will not be informative.
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:
# mod <- model.matrix(~0+samples$condition+samples$batch)
mod <- model.matrix(~0+samples$condition)

# make model terms easier to work with
colnames(mod) <- levels(samples$condition)

fit <- lmFit(log_counts, design=mod)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(samples$condition), 2))

comparisons <- list()
for(i in 1:nrow(condition_pairs)){
comparisons[[i]] <- as.character(condition_pairs[i,])
}

# vector to store de genes
sig_genes <- c()


# iterate over the contrasts, and perform a differential expression test for
# each pair
for (conds in comparisons) {
    # generate string contrast formula, "infLM24 - infLM4"
    contrast_formula <- paste(conds, collapse=' - ')

    contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
    contrast_fit <- contrasts.fit(fit, contrast_mat)
    eb <- eBayes(contrast_fit)

    # Grab highly ranked genes; this is a pretty stringent p-value cutoff, but
    # it serves to limit the total number of genes we will use for this
    # tutorial
    sig_genes <- union(sig_genes,
                       rownames(topTable(eb, number=Inf, p.value=0.05)))
}


# Filter out genes which were not differentially expressed for any contrast
log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]

log_counts_t <- t(log_counts)

# Construct similarity matrix
# determine the soft-thresholding power to use
powers <- c(c(1:10), seq(from = 14, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(log_counts_t, powerVector = powers, verbose = 5)
# Plot the results:

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
#png("4. soft-thresholding_power.png")

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.82,col="red")
#dev.off()


#png("5. mean_connectivity.png")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()


# construct adjacency matrix
softpower <- 14
adjacency.matrix <- adjacency(log_counts_t, power=softpower,
                                             type = "signed", corFnc = "cor")

# Turn the adjacency matrix to topologicaal overlap matrix to minimize
# the effects of noise and spurious associations
TOM <- TOMsimilarity(adjacency.matrix, TOMType = "signed")
dissTOM <- 1 - TOM


#set diagonal to NA to remove uninformative correlations
diag(adjacency.matrix) <- NA

# Adjacency matrix heatmap plot / network heatmap of selected genes
heatmap_indices <- sample(nrow(adjacency.matrix), 500) # sub-sample for visualization purposes

#png("6. adjacency_matrix_heatmap.png")
heatmap.2(t(adjacency.matrix[heatmap_indices, heatmap_indices]),
            col=redgreen(75),
            labRow=NA, labCol=NA, 
            trace='none', dendrogram='row',
            xlab='Gene', ylab='Gene',
            main='Adjacency matrix',
            density.info='none', revC=TRUE)
#dev.off()


# remove adjacency matrix and TOM to free up memory
rm(adjacency.matrix)
gc()


################################################################
## Detecting co-expression modules in R
################################################################

# view the dendrogram based on hierachical clustering of genes
gene_tree <- flashClust(as.dist(dissTOM), method = "average")

# plot the gene tree
#png("7. gene_tree.png")
#sizeGrWindow(12,9) #open graphical window
plot(gene_tree, xlab="", sub="", main = "Gene clustering based on TOM dissimilarity", 
     labels = FALSE, hang = 0.04)
#dev.off()

# identify the modules
module_labels <- cutreeDynamicTree(gene_tree, deepSplit = TRUE, 
                                   minModuleSize = 30)

#view
table(module_labels)

# convert labels to colours
module_colours <- labels2colors(module_labels)

# view
table(module_colours)

#### a list of 28 modules#
table(module_colours)

# visualize the gene tree and TOM matrix together using TOM plot
# if necessary, raise dissTOM to a power to make moderately strong connection more visible in heatmap
diag(dissTOM) <- NA

#png("8. gene_tree_and_dissTOM.png")
TOMplot(dissTOM, gene_tree, as.character(module_colours))
#dev.off()
# remove matrix to free memory
rm(dissTOM)
gc()


# plot gene dendrogram
#png(filename = "9. gene_tree_and_colours.png")
#sizeGrWindow(8,6) #open graphical window
plotDendroAndColors(gene_tree, module_colours, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colours")
#dev.off()

# get hub genes
# choose power 4: https://support.bioconductor.org/p/46342/
module_hub_genes <- chooseTopHubInEachModule(log_counts_t, module_colours, 
					     power = 4,type = "signed")
module_hub_genes
length(module_hub_genes)
##### A list of module hub genes#
module_hub_genes
        black          blue         brown          cyan     darkgreen
      "INTS6"       "PRKDC"       "MAPK1"      "DUOXA1"   "LOC788703"
     darkgrey    darkorange       darkred darkturquoise         green
      "RAB1A"       "KMT2E"       "LUC7L"        "ATRN"      "HNRNPM"
  greenyellow        grey60     lightcyan    lightgreen   lightyellow
      "DDX53"       "LHCGR"      "CHCHD3"       "SGMS1"    "ARHGAP35"
      magenta  midnightblue        orange          pink        purple
      "LRIG2"        "GJA5"       "ELANE"      "KANSL1"        "MPP6"
          red     royalblue        salmon       skyblue           tan
       "ABI1"       "CNOT4"     "MFSD14B"       "WDR33"       "TKDP1"
    turquoise         white        yellow
       "RHOF"        "ELF1"       "RC3H1"


#The section below is included for further checks so it may not be necessary to carry out this analysis.

# --------------------------------------------------------------------------------------------
# merge modules with very similar expression profiles as their genes are highly co-expressed
# get the module eigengenes
#module_eigengenes <- moduleEigengenes(log_counts_t, colors = module_colours)$eigengenes
#module_eigengenes
# calculate dissimilarity of module eigengenes using correlations
#module_eigengenes_diss <- 1 - cor(module_eigengenes)

# cluster module eigengenes
#module_eigengenes_tree <- flashClust(as.dist(module_eigengenes_diss), method = "average")

# choose height at which to cut the tree for merge i.e. the threshold
#module_eigengenes_thresh <- 0.25

# create plots for the results
#png("10. module_eigengenes_cluster.png")
#sizeGrWindow(7, 6)
#plot(module_eigengenes_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#abline(h=module_eigengenes_thresh, col="red")

#dev.off()

# merge the modules
#module_eigengenes_merge <- mergeCloseModules(log_counts_t, module_colours, 
 #                                            cutHeight = module_eigengenes_thresh)

# merged module colours
#merged_module_colours <- module_eigengenes_merge$colors
#merged_module_colours
#length(merged_module_colours)

# view
#table(merged_module_colours)

### a list of 13 modules


# eigengenes of new merged modules
#merged_module_eigengenes <- module_eigengenes_merge$newMEs


# plot the dendrogram with original and merged colours underneath
#sizeGrWindow(12, 9)
#png("11. merged-original_colours-original_dendro.png")
#plotDendroAndColors(gene_tree, cbind(module_colours, merged_module_colours), 
 #                   c("Dynamic Tree Cut", "Merged dynamic"), 
  #                  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
#dev.off()

# plot heatmap of eigengenes (orginal before merge)
#png("12. eigengenes_heatmap.png")
#plotEigengeneNetworks(module_eigengenes, "Eigengenes heatmap", marHeatmap = c(3,4,2,2),
 #                     plotDendrograms = FALSE, xLabelsAngle = 90)
#dev.off()



#-----------------------------------------------------------------------------------------
# rename some variables based on the module eigengene analysis for later use
#
# module colours
#module.colours <- merged.module.colours

# construct numerical labels corresponding to the colours
#colorOrder <- c("grey", standardColors(50))
#module.labels <- match(module.colours, colorOrder)-1

# module eigengenes
#module.eigengenes <- merged.module.eigengenes

# get hub genes
#merged_module_hub_genes <- chooseTopHubInEachModule(log_counts_t,
 #                                                   merged_module_colours,
  #                                                  power = 4,
   #                                                 type = "signed")

#### a list of merged module hub genes#
#merged_module_hub_genes 

##############################################################################
## Network export to cytoscape
##############################################################################

# select modules of interest

interesting_modules <- c('black', 'blue', 'brown','cyan',
                         'darkgreen','darkgrey','darkorange','darkred', 'darkturquoise','green',
                         'greenyellow','grey','grey60','lightcyan','lightgreen','lightyellow',
                         'magenta','midnightblue','orange','pink','purple','red','royalblue',
                         'salmon','skyblue','tan','turquoise','white','yellow') # all module colours, thus the whole network

# enriched modules
#enriched_modules <- c("black","tan","brown","blue","turquoise","magenta","darkturquoise",
 #                    "green","red","pink","salmon","lightyellow","purple","greenyellow")


# enriched modules
enriched_modules <- c("black","blue","brown","green","pink","red","turquoise","yellow")


# obtain gene ids
gene_ids <- rownames(log_counts)

# select module genes
#inModules <- is.finite(match(module.colours, interesting.modules)) # whole network modules

#inModules <- is.finite(match(module.colours, enriched.modules)) # enriched modules
inModules <- is.finite(match(module_colours, c("black", "blue", "brown", "green", "pink", "red", "turquoise", "yellow"))) # individual modules

modGenes <- gene_ids[inModules]

# select the corresponding dissTOM based on module genes
modTOM <- TOM[inModules, inModules]
dimnames(modTOM) <- list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
exportNetworkToCytoscape(modTOM,
                         edgeFile = "CytoscapeInput-edges_skyblue_white_module_thresh0.txt",
                         nodeFile = "CytoscapeInput-nodes_skyblue_white_module_thresh0.txt",
                         weighted = TRUE,
                         threshold = 0,
                         nodeNames = modGenes,
                         nodeAttr = module_colours[inModules]);


# Also, export the network as graphml format
# use export_network_to_graphml function

source("network_export_graphml.R")

# the whole network
entire_network <- export_network_to_graphml(TOM, 
                                            filename = "entire_network_thresh0.graphml",
                                            threshold = 0.4, #### why not 0.4
                                            nodeAttr = gene_ids,
                                            max_edge_ratio = 3)
                                            
 
# network modules
# create a dataframe with node attributes
interesting_module_colours <- module_colours[inModules] #get enriched module colours from module.colours
node_attributes <- cbind(modGenes, module=interesting_module_colours) # node atrr. for enriched modules

#node.attributes <- cbind(modGenes, module=module.colours) # get node attr. for whole network
node_attributes <- as.data.frame(node_attributes)

# Add RGB versions of colour modules
node_attributes$colourRGB <- col2hex(node_attributes$module)

modules_network <- export_network_to_graphml(modTOM, 
                                             filename = "modules_network_thresh0.02.graphml",
                                             threshold = 0.02, ###threshold change ######
                                             nodeAttrDataFrame = node_attributes)

# write out a node attributes files with hexadecimal colour names for module genes
write.table(node_attributes, 
            file = "Cytoscape_node_attributes_interesting_modules_skyblue_white.txt",
            row.names = FALSE, 
            quote = FALSE, sep = "\t")

## Functional Analysis

### Loading annotations from the packages





