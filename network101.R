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

# Make sure results are reproducible
set.seed(1)

# Read sample meta-data file
samples <- read.csv('/home/user/Desktop/SraRunTable.csv', header=TRUE)
samples<-samples[,c(1,16,25)] # ensures only 3 columns are left

#see results after retaining most important columns (sample name and serologic_response_status)
kable(samples)


# Read in count data
raw_counts <- read.table('/home/user/Desktop/NA.txt', row.names=1, header=TRUE)
head(raw_counts)
dim(raw_counts)

names(raw_counts)
raw_counts<-raw_counts[-c(1:5),]


# Data Preparation
# Sample check
# add a colorbar along the heatmap with sample condition

num_conditions <- nlevels(samples$condition)
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(samples$condition)]

heatmap.2(cor(raw_counts), RowSideColors=cond_colors,
          trace='none', main='Raw Correlations')

#Low count filtering
# Remove all rows with less than 1 count per sample, at least 45 reads per gene
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)

sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))

# [1] "Removing 11834 low-count genes (17087 remaining)."

# Log2 transformation
log_counts <- log2(raw_counts + 1)

# visualize the data as density plots, like histogram for each sample

x = melt(as.matrix(log_counts))

colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + geom_density()

heatmap.2(cor(log_counts), RowSideColors=cond_colors,
          trace='none', main='log2-transformed correlations)')

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




#########################################################################################################################################################################################################################
#                                                  Co-expression network construction
#####################################################################################################################################################################
# Construct similarity matrix

#'
#' Similarity measure which combines elements from Pearson correlation and
#' Euclidean distance.
#' 
cordist <- function(dat) {
    cor_matrix  <- cor(t(dat))

    dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
    dist_matrix <- log1p(dist_matrix)
    dist_matrix <- 1 - (dist_matrix / max(dist_matrix))

    sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}


sim_matrix <- cordist(log_counts)

heatmap_indices <- sample(nrow(sim_matrix), 500)

heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
            col=redgreen(75),
            labRow=NA, labCol=NA, 
            trace='none', dendrogram='row',
            xlab='Gene', ylab='Gene',
            main='Similarity matrix',
            density.info='none', revC=TRUE)

# to construct an adjacency matrix, the similarity matrix needs to be transposed 

log_counts_t <- t(log_counts)


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

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()


#####################################################################################################################

##################################################################################################################


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

