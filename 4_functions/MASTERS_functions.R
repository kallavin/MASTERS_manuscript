# WGCNA-Specific

getSize <- function(ExprSet, Network, module)
{
  modLabels = Network$colors
  modColors = labels2colors(Network$colors)
  ME <- names(ExprSet)[modColors==module]
  length(ME)
}

NetworkExam <- function(ExprSet, modColors){
  MEs1 = moduleEigengenes(ExprSet, modColors)$eigengenes
  MEs = orderMEs(MEs1)
  
  distance = 1-(1+cor(MEs,use="p"))/2
  MEcluster = hclust(as.dist(distance), method = "average")
  
  plot(MEcluster, main = "Sample clustering of Module Eigengenes", sub="", xlab="", cex.lab = 1,cex.axis = 1, cex.main = 1)
  
  MDS = cmdscale(as.dist(distance),2) 
  plot(MDS, col=modColors)
}

ModuleTable <- function(ExprSet, Network, modColors){
  MEs = moduleEigengenes(ExprSet, modColors)$eigengenes
  MEs = orderMEs(MEs)
  modNames <- substring(names(MEs),3) #drop ME
  
  geneModuleMembership = as.data.frame(cor(ExprSet, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(ExprSet)));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  #Make a table with the columns designating the module name
  modsTable <- data.frame(row.names= modNames)
  modsTable$size <- NA
  modsTable$kME <- NA
  
  for (i in 1:length(MEs)) {
    module=modNames[i]
    modsTable$size[i] <- getSize(ExprSet = ExprSet,Network = Network, module = module)
  }
  
  for (i in 1:length(MEs)){
    module=modNames[i]
    column = match(module, modNames);
    moduleGenes = modColors==module;
    modsTable$kME[i] <-round(mean(geneModuleMembership[moduleGenes, column]),4)
  }
  modsTable
}

createTablefor <- function(module, ExprSet, modColors, trait, filepath){
  trait = as.data.frame(trait)
  MEs = moduleEigengenes(ExprSet, modColors)$eigengenes
  modNames = substring(names(MEs), 3)
  nSamples = nrow(ExprSet)
  
  geneModuleMembership = as.data.frame(cor(ExprSet, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(ExprSet, trait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(trait), sep="");
  names(GSPvalue) = paste("p.GS.", names(trait), sep="");
  
  column = match(module, modNames);
  moduleGenes = modColors==module;
  
  modmem <- data.frame(abs(geneModuleMembership[moduleGenes, column]),row.names = rownames(geneModuleMembership[moduleGenes,]))
  
  names(modmem)[names(modmem)=="abs.geneModuleMembership.moduleGenes..column.."] <- "membership"
  
  modmem <- modmem %>%
    #mutate("traitsig" = abs(geneTraitSignificance[moduleGenes, 1]), 
    mutate("geneIDs" = row.names(modmem)) %>% 
    separate(geneIDs, into = c("ensg_version","symbol"), sep="_")
  
  row.names(modmem) <- modmem$geneIDs
  modmem$geneIDs <- NULL
  filename = paste(c(module,"GeneTable.csv"), collapse="")
  write.csv(modmem, file=paste0(filepath, filename))
}

createPlotfor <- function(module, ExprSet, modColors, Network, trait, filepath){
  trait = as.data.frame(trait)
  MEs = moduleEigengenes(ExprSet, modColors)$eigengenes
  modNames = substring(names(MEs), 3)
  nSamples = nrow(ExprSet)
  
  #calculate most related genes within module
  geneModuleMembership = as.data.frame(cor(ExprSet, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  geneModuleMembership$geneID <- rownames(geneModuleMembership)
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  #calculate genes with strongest association to trait
  geneTraitCor = as.data.frame(cor(ExprSet, trait, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor$trait), nSamples));
  geneTraitCor$geneID <- rownames(geneTraitCor)
  names(GSPvalue) = paste("p.GS.", names(trait), sep="");
  names(geneTraitCor) = c(paste("gene_cor_", names(trait), sep=""), "geneID");
  
  #identify where in the full module list the genes in the module of interest fall
  column = match(module, modNames);
  moduleGenes = modColors==module;
  
  #prepare to join two data frames
  modmem <- data.frame(geneModuleMembership[moduleGenes, column],row.names = rownames(geneModuleMembership[moduleGenes,]))
  names(modmem)[names(modmem)=="geneModuleMembership.moduleGenes..column."] <- "membership"
  modmem$geneID <- rownames(modmem)
  
  modmem <- left_join(modmem, geneTraitCor, by="geneID") 
  modmem <- modmem %>%
    dplyr::select("geneID","membership","gene_cor_trait")
  
  #establish criteria for cutoff into the fractions of interest
  medianmembership <- median(modmem$membership, na.rm = TRUE)
  mediantraitsig <- median(modmem$gene_cor_trait, na.rm = TRUE)
  
  #need to identify which transcripts have modmem that fall into the top ten percent in terms of module membership
  highmodmem <- modmem[order(-modmem$membership),]
  tenpercentMM <- round(nrow(modmem[-rank(modmem$membership),])/10)
  
  top_half <- modmem[modmem$membership >= medianmembership,]
  quadrant4 <- modmem[(modmem$membership >= medianmembership) & (modmem$trait >= abs(mediantraitsig)),] #use abs() since sometimes the correlation will be opposite with trait
  highkME <- modmem[modmem$membership >= 0.8,]
  
  modmem <- modmem %>%
    mutate("top_half"= NA,"quadrant_4" = NA,"kME_above_0.8"= NA,"top_ten_percent"= NA, "top_ten_percent_highGS"= NA, "gene_biotype"=NA) %>%
    filter(!is.na(gene_cor_trait))
  
  #define which fractions each transcript is in 
  for (i in 1:length(rownames(modmem))){
    if (modmem$membership[i] >= medianmembership){
      modmem$top_half[i]="TRUE"
      if (modmem$gene_cor_trait[i] >= abs(mediantraitsig)){
        modmem$quadrant_4[i]="TRUE"
      }
    }
    if (modmem$membership[i] >= 0.8){
      modmem$kME_above_0.8[i] = "TRUE"
    }
    if (round(rank(-modmem$membership))[i] <= round(nrow(modmem)/10)){
      modmem$top_ten_percent[i]="TRUE"
      if (modmem$gene_cor_trait[i] >= abs(mediantraitsig)){
        modmem$top_ten_percent_highGS[i]="TRUE"
      }
    }
    
  }
  
  #annotate the transcripts with function and name if available using annotmods
  load(file="/data/project/bamman-lab/Parkinsons/ExercisePLIER/2_pipeline/annotationFile.csv")
  for (i in 1:nrow(modmem)){
    ID <-modmem$geneID[i]
    modmem$gene_biotype[i] <- annotationFile$gene_biotype[match(ID, annotationFile$external_gene_name)]
  }
  
  filename <- paste(module,"_module_details.csv", sep="")
  write.csv(modmem, file = paste("/data/project/bamman-lab/MASTERS/3_output/1_MASTERS_WGCNA/GeneDetails/paired", filename))
  
  sizeGrWindow(7,7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                     abs(geneTraitCor[moduleGenes, 1]),
                     xlab = paste("Gene Membership in", module, "Module"), 
                     ylab = paste("Gene correlation for Trait"), 
                     cex.main = 1, cex.lab = .7, cex.axis = 1, bg = module, col = module, pch = 21)
  filename <- paste("Scatterplot", module,".pdf", sep="", 
                    abline(h=abs(mediantraitsig), col="black"), 
                    abline(v=abs(medianmembership), col="black"))
  
  dev.copy2pdf(file=paste0(filepath, filename), width=6, height=4)
} 


#I will create a function to output the information necessary for visualization in Cytoscape. 
getCytoscape <- function(module, detailTable, ExprSet, fraction, modColors, TOM_matrix){
  #user must always input the module and the fraction with quotes!
  
  #identify the fraction of interest based on user input
  if (fraction =="top_half"){
    detailTable <- detailTable[!is.na(detailTable$top_half),]
  }
  if (fraction =="quadrant_4"){
    detailTable <- detailTable[!is.na(detailTable$quadrant_4),]
  }
  if (fraction=="kME_above_0.8"){
    detailTable <- detailTable[!is.na(detailTable$kME_above_0.8),]
  }
  if (fraction=="top_ten_percent"){
    detailTable <- detailTable[!is.na(detailTable$top_ten_percent),]
  }
  if (fraction=="top_ten_percent_highGS"){
    detailTable <- detailTable[!is.na(detailTable$top_ten_percent_highGS),]
  }
  if (fraction=="all"){
    detailTable <- detailTable
  }
  
  #subset to gene names for the protein-coding genes only
  subset <- detailTable %>%
    separate(geneID, into = c("ensg_version","ID"), sep = "_",remove = FALSE) %>%
    dplyr::select(geneID)
  
  #subset the large matrix to define what meets criteria
  all_genes <- colnames(ExprSet) 
  inSubset = is.finite(match(all_genes, subset$geneID))
  subset_ids = all_genes[inSubset];
  subset_ids_gene <- gsub(pattern = ".*_","",subset_ids)
  
  #which enables subsetting the TOM matrix to only what meets criteria
  subset_TOM = TOM_matrix[inSubset, inSubset];
  dimnames(subset_TOM) = list(subset_ids_gene, subset_ids_gene)
  
  #perform export to Cytoscape
  cyt = exportNetworkToCytoscape(subset_TOM,
                                 edgeFile = paste("../../3_output/1_MASTERS_WGCNA/Cytoscape/GeneBiotype-edges-", 
                                                  paste(module, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("../../3_output/1_MASTERS_WGCNA/Cytoscape/GeneBiotype-nodes", 
                                                  paste(module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = subset_ids_gene,
                                 altNodeNames = subset_ids_gene,
                                 nodeAttr = modColors[inSubset]);
}


# PLIER-Specific
##Chikina tscale code
tscale=function(mat){
  t(scale(t(mat)))
}

##plot the heatmaps
customplotTopZcontinuous <- function (plierRes, variable, priorMat, top = 50, index = NULL, regress = F, allLVs = F, ...) 
{
  continuousdf <- continuousdf[order(continuousdf[,variable]),]
  plotdf <- continuousdf %>% 
    dplyr::select("Batch",variable) 
  
  data <- yN[,c(rownames(continuousdf))]
  
  target <- rownames(plotdf)
  data <- data[,match(target, colnames(data))]
  
  #subset data to basal only for this heatmap
  data = data[rownames(plierRes$Z), ] #genes
  
  #default code from MSSM
  priorMat = priorMat[rownames(plierRes$Z), ]
  ii = index
  tmp = apply(-plierRes$Z[, ii, drop = F], 2, rank)
  nn = character()
  nncol = character()
  nnpath = character()
  nnindex = double()
  
  for (i in 1:length(ii)) {
    nn = c(nn, nntmp <- names(which(tmp[, i] <= top)))
    nncol = c(nncol, rep(rownames(plierRes$B)[ii[i]], length(nntmp)))
    nnpath = c(nnpath, rowSums(priorMat[nntmp, plierRes$U[, 
                                                          ii[i]] > 0, drop = F]) > 0)
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  names(nncol) = nn
  nncol = strtrim(nncol, 30)
  nnrep = names(which(table(nn) > 1))
  
  if (length(nnrep) > 0) {
    nnrep.im = match(nnrep, nn)
    nn = nn[-nnrep.im]
    nncol = nncol[-nnrep.im]
    nnpath = nnpath[-nnrep.im]
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  nnpath[nnpath == "TRUE"] = "inPathway"
  nnpath[nnpath == "FALSE"] = "notInPathway"
  nncol = as.data.frame(list(nncol, nnpath))
  names(nncol) = c("pathway", "present")
  
  toPlot = tscale(data[nn, ])
  if (regress) {
    for (i in ii) {
      gi = which(nnindex == i)
      toPlot[gi, ] = resid(toPlot[gi, ], model.matrix(~t(plierRes$B[-i, 
                                                                    colnames(toPlot)])))
    }
  }
  maxval = max(abs(toPlot))
  
  #adjust for data set
  BatchCol <- c("darkgreen","deepskyblue","tan1")
  plotdf$Batch <- factor(plotdf$Batch, levels = c("1","2","3"))
  names(BatchCol) <- levels(plotdf$Batch)
  
  Variable <- names(plotdf[names(plotdf)==variable])
  pretty <- brewer.pal(9, "YlOrRd")
  newcol <- colorRampPalette(pretty)
  ContinCols <- newcol(100)
  
  ll = c(inPathway = "black", notInPathway = "beige")
  
  AnnColour <- list(Batch = BatchCol, Variable = ContinCols, present = ll)
  AnnColour[[variable]]<- AnnColour[["Variable"]]
  
  pheatmap(toPlot, breaks = seq(-maxval, maxval, length.out = 99), 
           color = colorpanel(100, "blue", "white", "red"), annotation_row = nncol, 
           show_colnames = F, annotation_colors = AnnColour, cluster_cols = F, cellwidth = 7, cellheight= 10, 
           annotation_col = plotdf, ...)
}

customplotTopZ <- function(plierRes, data, priorMat, top = 10, index = NULL, regress = F, allLVs = F, ...) 
{
  #this subsets the data from newyN to only the data that were handled by PLIER, i.e. going from 17747 to 5884
  data = data[rownames(plierRes$Z), ]
  
  #this subsets the pathway data to only the pathways PLIER found genes in 
  priorMat = priorMat[rownames(plierRes$Z), ]
  
  #ultimately this sets the index to the LV of interest, since we routinely will have allLVs set to TRUE
  ii = which(colSums(plierRes$U) > 0)
  if (!allLVs) {
    if (!is.null(index)) {
      ii = intersect(ii, index)
    }
  }
  else {
    ii = index
  }
  
  #here we are taking the negative of gene contributions to LVs and apply the rank functions to columns
  #it is probably made negative to keep the largest values the top ranks
  #the drop argument helps us keep the correct dimensions of the plierResult$Z
  #the rank function returns the ranks, i.e., where the values fall in relation to the others in the vector
  tmp = apply(-plierRes$Z[, ii, drop = F], 2, rank)
  
  #initializing vectors to hold the information (double is a numeric non-integer vector)
  nn = character()
  nncol = character()
  nnpath = character()
  nnindex = double()
  
  #the for loop below first identifies gene names from the tmp that are in the top n (user-defined category)
  #next it looks through plierResult$B and gets LV associations for genes and LV pathway association if one was found
  #next it will look through the allPaths and get the pathway associations for those genes (T/F)
  #finally it populates a vector defining the LV association for all the genes to plot 
  for (i in 1:length(ii)) {
    nn = c(nn, nntmp <- names(which(tmp[, i] <= top)))
    nncol = c(nncol, rep(rownames(plierRes$B)[ii[i]], length(nntmp)))
    nnpath = c(nnpath, rowSums(priorMat[nntmp, plierRes$U[, 
                                                          ii[i]] > 0, drop = F]) > 0)
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  #here the 50 genes are placed into a vector, and the pathway/LV associations are trimmed to a length of 30
  #it will scan to see if any genes are replicated in the data set.
  names(nncol) = nn
  nncol = strtrim(nncol, 30)
  nnrep = names(which(table(nn) > 1))
  
  #if any duplicates were found, it drops them
  if (length(nnrep) > 0) {
    nnrep.im = match(nnrep, nn)
    nn = nn[-nnrep.im]
    nncol = nncol[-nnrep.im]
    nnpath = nnpath[-nnrep.im]
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  #define desired colors (amended by KL)
  TimepointCol <- c("white","grey","black")
  simpledf$TimePoint <- factor(simpledf$TimePoint, levels = c("1","2","3"))
  names(TimepointCol) <- levels(simpledf$TimePoint)
  
  BatchCol <- c("darkgreen","deepskyblue","tan1")
  simpledf$Batch <- factor(simpledf$Batch, levels = c("1","2","3"))
  names(BatchCol) <- levels(simpledf$Batch)
  
  SexCol <- c("purple","turquoise2")
  simpledf$Sex <- factor(simpledf$Sex, levels = c("1","2"))
  names(SexCol) <- levels(simpledf$Sex)
  
  ll = c(inPathway = "green", notInPathway = "beige")
  AnnColour <- list(TimePoint = TimepointCol, Sex = SexCol, Batch = BatchCol, present = ll)
  
  nnpath[nnpath == "TRUE"] = "inPathway"
  nnpath[nnpath == "FALSE"] = "notInPathway"
  
  #the LV associations and pathway associations vectors are merged to allow plotting 
  nncol = as.data.frame(list(nncol, nnpath))
  names(nncol) = c("pathway", "present")
  
  #recall tscale scales a transposed matrix, in this case the normalized expression for the top LV genes
  #the scale function centers the matrix, then t transposes it back
  toPlot = tscale(data[nn, ])
  if (regress) {
    for (i in ii) {
      gi = which(nnindex == i)
      toPlot[gi, ] = resid(toPlot[gi, ], model.matrix(~t(plierRes$B[-i, 
                                                                    colnames(toPlot)])))
    }
  }
  #sets the colors to be specific to the values in THAT particular plot
  maxval = max(abs(toPlot))
  
  #modified pheatmap to allow gaps, annotation pulling from simpledf, and desired uniform margins
  pheatmap(toPlot, breaks = seq(-maxval, maxval, length.out = 99), 
           color = colorpanel(100, "blue", "white", "red"), annotation_row = nncol, 
           show_colnames = F, annotation_colors = AnnColour, cluster_cols = F, cellwidth = 7, cellheight= 10,...)
}

##Boxplots
PlierBox <- function(plierRes, index = NULL, subject.data, group.var, paired, subject.var, ...)
{
  #user defines plierResult data and LV of interest (index), subject data, and how to group them
  #function pulls data from plierResult$B
  score <- data.frame(row.names = colnames(plierRes$B)) %>% 
    rownames_to_column(var = "seqID") %>% 
    mutate(avgExpr = plierRes$B[index,])
  
  BoxPlotter <- subject.data %>% 
    mutate(avgExpr = score$avgExpr[match(seqID, score$seqID)])
  
  #this part can all stay the same unless the user defines paired/non paired
  gg.group <- subject.data[,group.var]
  
  limits <- c(min(BoxPlotter$avgExpr)-.1, max(BoxPlotter$avgExpr)+.1)
  
  if (paired){
    q <- ggplot(BoxPlotter) + 
      geom_boxplot(aes(x = BoxPlotter[,group.var], y = avgExpr), 
                   outlier.shape = NA, color="black", fill=c("white","darkgrey"), lwd = .75) + 
      geom_point(aes(x = BoxPlotter[,group.var], y = avgExpr, group = BoxPlotter[,subject.var]),
                 size=2, position = position_dodge(width = 0.1)) + 
      geom_line(aes(x = BoxPlotter[,group.var], y = avgExpr, group = BoxPlotter[,subject.var]), 
                position = position_dodge(width = 0.1), size = .75) + 
      scale_y_continuous(name = paste(c("Values for LV", index), collapse=""), limits=limits, 
                         breaks = seq(round(limits[1],1), round(limits[2],1),1)) +
      theme_classic()+  
      theme(plot.title = element_blank(), 
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 20, colour = "black", face = "bold"), 
            axis.text.y = element_text(size = 16, colour = "black", face = "bold"), 
            axis.text.x = element_text(size = 20, colour = "black", face = "bold"), 
            axis.ticks.y = element_line(size = 1, colour = "black"),
            axis.ticks.x = element_line(size = 1, colour = "black"),
            axis.ticks.length = unit(5, "pt"),
            axis.line = element_line(colour = 'black', size = 1))
  }
  
  else if (!paired){
    q <- ggplot(BoxPlotter) + 
      geom_boxplot(aes(x = BoxPlotter[,group.var], y = avgExpr), 
                   outlier.shape = NA, color="black", fill=c("white","darkgrey"), lwd = .75) + 
      geom_point(aes(x = BoxPlotter[,group.var], y = avgExpr, group = BoxPlotter[,subject.var]),
                 size=2, position = position_dodge(width = 0.1)) + 
      
      scale_y_continuous(name = paste(c("Values for LV", index), collapse=""), limits=limits, 
                         breaks = seq(round(limits[1],1), round(limits[2],1),1)) +
      theme_classic()+  
      theme(plot.title = element_blank(), 
            axis.title.x.bottom = element_blank(),
            axis.title.y = element_text(size = 20, colour = "black", face = "bold"), 
            axis.title.y = element_blank(), 
            axis.text.y = element_text(size = 16, colour = "black", face = "bold"), 
            axis.text.x = element_text(size = 20, colour = "black", face = "bold"), 
            axis.ticks.y = element_line(size = 1, colour = "black"),
            axis.ticks.x = element_line(size = 1, colour = "black"),
            axis.ticks.length = unit(5, "pt"),
            axis.line = element_line(colour = 'black', size = 1))
  }
  
}
## Tables
LVtablemaker <- function(plierRes, data, priorMat, index=NULL, top=NULL, path=NULL,...)
{
  data = data[rownames(plierRes$Z), ]
  priorMat = priorMat[rownames(plierRes$Z), ]
  ii = index
  tmp = apply(-plierRes$Z[, ii, drop = F], 2, rank)
  
  nn = character()
  nncol = character()
  nnpath = character()
  nnindex = double()
  
  for (i in 1:length(ii)) {
    nn = c(nn, nntmp <- names(which(tmp[, i] <= top)))
    nncol = c(nncol, rep(rownames(plierRes$B)[ii[i]], length(nntmp)))
    nnpath = c(nnpath, rowSums(priorMat[nntmp, plierRes$U[, 
                                                          ii[i]] > 0, drop = F]) > 0)
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  names(nncol) = nn
  nncol = strtrim(nncol, 30)
  nnrep = names(which(table(nn) > 1))
  
  if (length(nnrep) > 0) {
    nnrep.im = match(nnrep, nn)
    nn = nn[-nnrep.im]
    nncol = nncol[-nnrep.im]
    nnpath = nnpath[-nnrep.im]
    nnindex = c(nnindex, rep(ii[i], length(nntmp)))
  }
  
  nnpath[nnpath == "TRUE"] = "inPathway"
  nnpath[nnpath == "FALSE"] = "notInPathway"
  
  nncol = as.data.frame(list(nncol, nnpath))
  names(nncol) = c("pathway", "present")
  
  #creates the table using toPlot
  
  LVtable <- nncol
  temptable <- getBM(attributes=c('external_gene_name','description','gene_biotype'), 
                     filters = 'external_gene_name', values = rownames(LVtable), mart = ensembl)
  LVtable$description <- temptable$description[match(rownames(LVtable),temptable$external_gene_name)]
  LVtable$biotype <- temptable$gene_biotype[match(rownames(LVtable),temptable$external_gene_name)]
  
  #drop bracketed bit, using forward slash to escape the regular expression syntax
  LVtable$description <- as.character(LVtable$description)
  LVtable$description <- gsub("\\[.*","", LVtable$description)
  
  filename <- paste(c(path,"/LV",index,"_table_Top",top,".csv"), collapse="")
  
  head(LVtable)
  write.csv(LVtable, file=filename)
}


# Integration and Other
mod_LV_check = function(plierResult, genetable,... ){
  z <- plierResult$Z
  colnames(z) <- rownames(plierResult$B)
  
  genetable <- genetable %>% 
    filter(symbol %in% rownames(plierResult$Z)) #108 out of 307
  
  filt_z <- as.data.frame(z) %>% #hopefully more computationally efficient
    filter(rownames(z) %in% genetable$symbol)
  try <- t(filt_z)
  
  mod_LV_associations <- as.data.frame(genetable$symbol) %>% 
    rename("geneID" = "genetable$symbol") %>% 
    mutate("LV_number" = "LV") %>% 
    mutate("LV_path" = "LVpath") %>% 
    arrange(match(geneID,rownames(filt_z)))
  
  for (i in 1:nrow(mod_LV_associations)){
    mod_LV_associations$LV_number[i] = which.max(filt_z[i,])
    mod_LV_associations$LV_path[i] = names(filt_z)[which.max(filt_z[i,])]
  }
  
  mod_LV_associations
}

#a function to look for matches in given lists
IDoverlaps <- function(genelist1, genelist2){
  inList <- is.finite(match(genelist1, genelist2))
  matches <- genelist1[inList]
  
  if ((length(matches)) >= 0){
    matches
  } else {
    print("No matches found.")
  }
}

