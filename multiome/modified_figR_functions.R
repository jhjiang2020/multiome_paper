# Script for figR helper functions
## https://github.com/buenrostrolab/FigR
## modified to allow customized analysis

library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(networkD3)
library(GGally)
library(igraph)
library(network)
library(dplyr)


### modified figR code
smoothScoresNN <- function(NNmat,
                           TSSmat,
                           geneList = NULL,
                           barcodesList=NULL,
                           nCores = 1)
{
  if (is.null(rownames(NNmat)))
    stop("NN matrix has to have matching cell IDs as rownames\n")
  if (!all.equal(rownames(NNmat), colnames(TSSmat)))
    stop("Nearest-neighbor matrix and TSS activity score matrix don't have matching cells ..\n")
  cat("Number of cells in supplied TSS matrix: ", ncol(TSSmat),
      "\n")
  cat("Number of genes in supplied TSS matrix: ", nrow(TSSmat),
      "\n")
  cat("Number of nearest neighbors being used per cell for smoothing: ",
      ncol(NNmat), "\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% rownames(TSSmat)))) {
      cat("One or more of the gene names supplied is not present in the TSS matrix provided: \n")
      cat(geneList[!geneList %in% rownames(TSSmat)], sep = ", ")
      cat("\n")
      stop()
    }
    cat("Running TSS score smoothing for genes:", geneList,
        sep = "\n")
    cat("........\n")
    TSSmat <- TSSmat[rownames(TSSmat) %in% geneList, ]
  }
  else {
    if(nrow(TSSmat) > 10000){
      cat("Running smoothing for all genes in TSS matrix! (n = ",
          nrow(TSSmat), ") This is bound to take more time than querying specific markers ..\n",
          sep = "")
    }
  }
  opts <- list()
  pb <- txtProgressBar(min = 0, max = ncol(TSSmat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  
  
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  
  if(!is.null(barcodesList)){
    cat("Subsetting to ",length(barcodesList)," barcodes in dataset..\n")
    NNmat <- NNmat[barcodesList,]
  }
  cat("Running in parallel using ", nCores, "cores ..\n")
  matL <- foreach::foreach(x=1:nrow(NNmat),.options.snow = opts,.packages = c("Matrix","data.table","dplyr")) %dopar% {
    smoothedScore <- data.table(Matrix::rowMeans(TSSmat[, NNmat[x,]]))
    rownames(smoothedScore) <- rownames(TSSmat)
    colnames(smoothedScore) <- rownames(NNmat)[x]
    smoothedScore
  }
  
  parallel::stopCluster(cl)
  
  close(pb)
  cat("Merging results ..\n")
  smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix(sparse=TRUE)
  rownames(smoothedMat) <- rownames(TSSmat)
  #stopifnot(all.equal(colnames(smoothedMat), colnames(TSSmat)))
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed),
            "\n"))
  
  return(smoothedMat)
}

runFigR.2 <- function(ATAC.se, # SE of scATAC peak counts. Needed for chromVAR bg peaks etc.
                      dorcK=30, # How many dorc kNNs are we using to pool peaks
                      dorcTab, # peak x DORC connections (should contain indices relative to peaks in ATAC.se)
                      n_bg=50, # No. of background peaks to use for motif Z test
                      genome, # One of mm10, hg19, hg38, with no default
                      dorcMat, # Expect smoothed
                      rnaMat, # Expect smoothed
                      dorcGenes=NULL, # If only running on a subset of genes
                      nCores=1
){
  source("~/Oxford/Code/R/multiome/004_Downstream_figr/utils.R", local = T)
  # Must be matched data
  stopifnot(all.equal(ncol(dorcMat),ncol(rnaMat)))
  
  # Expects "Gene" / "Peak" in dorcTab
  if(!all(c("peak","gene") %in% colnames(dorcTab)))
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")
  
  if(all(grepl("chr",dorcTab$peak,ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    
    if(!(all(grepl("chr",rownames(ATAC.se),ignore.case = TRUE))))
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")
    
    if(!all(dorcTab$peak %in% rownames(ATAC.se)))
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  } else{
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
    # If using index, make sure no indices are outside range of SE
    if(max(dorcTab$peak) > nrow(ATAC.se))
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }  
  
  
  if(is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  } else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], sep = ", ")
      cat("\n")
      stop()
    }
  }
  
  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
  rownames(DORC.knn) <- rownames(dorcMat)
  
  if (is.null(rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
  }
  
  if(grepl("hg",genome)){
    #pwm <- chromVARmotifs::human_pwms_v2
    pwm <- readRDS("~/Oxford/Code/R/multiome/004_Downstream_figr/cisBP_human_pfms_2021.rds")
  } else {
    #pwm <- BuenRTools::mouse_pwms_v3
    pwm <- readRDS("~/Oxford/Code/R/multiome/004_Downstream_figr/cisBP_mouse_pfms_2021.rds")
  }
  
  # Old motif naming convention
  if(all(grepl("_",names(pwm),fixed = TRUE)))
    names(pwm) <- extractTFNames(names(pwm))
  
  message("Removing genes with 0 expression across cells ..\n")
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat)!=0,]
  myGeneNames <- gsub(x = rownames(rnaMat),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
  rownames(rnaMat) <- myGeneNames
  
  # Only non-zero expression TFs (also found in rnaMat)
  motifsToKeep <- intersect(names(pwm),myGeneNames)
  
  # This has to be done on the full SE (same peakset used as input to dorc calling)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se,pwms = pwm[motifsToKeep],genome=genome)
  
  # Keep TFs with some peak x motif match
  motif_ix <- motif_ix[,Matrix::colSums(assay(motif_ix))!=0] 
  
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat,1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.mat),rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  } else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }
  
  # For each DORC, do motif enrichment among dorc sig Peaks, and correlation of DORC accessibility (smoothed) to TF RNA levels
  
  cat("Testing ",length(motifsToKeep)," TFs\n")
  cat("Testing ",nrow(dorcMat)," DORCs\n")
  library(doParallel)
  if(nCores > 1)
    message("Running FigR using ",nCores," cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g=dorcGenes,
                         .options.snow = opts, 
                         #"BuenRTools"
                         .packages = c( "dplyr","Matrix")) %dopar%   {
                           # Take peaks associated with this DORC gene and its k neighbors
                           # Pool and use union for motif enrichment
                           DORCNNpeaks <- unique(dorcTab$peak[dorcTab$gene %in% c(g,dorcGenes[DORC.knn[g,]])])
                           
                           if(usePeakNames)
                             DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks) # Convert to index relative to input
                           
                           mZ <- motifPeakZtest(peakSet = DORCNNpeaks,
                                                bgPeaks = bg,
                                                tfMat = assay(motif_ix))
                           mZ <- mZ[,c("gene","z_test")]
                           colnames(mZ)[1] <- "Motif"
                           colnames(mZ)[2] <- "Enrichment.Z"
                           mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
                           mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
                           mZ <- cbind("DORC"=g,mZ)
                           # Correlate smoothed dorc with smoothed expression, with spearman
                           corr.r <- cor(dorcMat[g,],t(as.matrix(rnaMat[mZ$Motif,])),method = "spearman")
                           stopifnot(all.equal(colnames(corr.r),mZ$Motif))
                           mZ$Corr <- corr.r[1,] # Correlation coefficient
                           mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
                           mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
                           mZ$Corr.log10P <- sign(mZ$Corr.Z)*-log10(mZ$Corr.P)
                           return(mZ)
                         }
  cat("Finished!\n")
  cat("Merging results ..\n")
  # Merge and save table for downstream filtering and plotting (network)
  TFenrich.d <- do.call('rbind',mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL
  
  # Make combined score based on multiplication
  # Here, we only sign by corr
  # Since sometimes we lose digit precision (1 - v small number is 1, instead of 0.9999999..)
  # Use Rmpfr, increase precision limits above default (100 here)
  TFenrich.d <- TFenrich.d %>% mutate("Score"=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(Enrichment.P,100))*(1-Rmpfr::mpfr(Corr.P,100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  TFenrich.d
}

# Bar plot of ranked activator / repressors
rankDrivers <- function(figR.d, pct=0.05, 
                        myLabels=NULL){
  figR.summ <- figR.d %>%  group_by(Motif) %>% 
    summarise(Score=mean(Score)) %>% 
    arrange(Score) %>% 
    mutate(Motif=factor(Motif,levels=as.character(Motif)))
  
  # Top and bottom %ile labels
  figR.summ$Label <- as.character(figR.summ$Motif)
  
  if(is.null(myLabels)){
    # Use quantiles to define what labels are drawn
    figR.summ$Label[figR.summ$Score >= quantile(figR.summ$Score,pct) & 
                      figR.summ$Score <= quantile(figR.summ$Score,1-pct)] <- ""
  } else {
    # Only highlight user-specified
    figR.summ$Label[!figR.summ$Label %in% myLabels] <- ""
  }
  
  
  gAll <- ggplot(figR.summ,aes(x=Motif,y=Score,label=Label)) + 
    geom_bar(size=0.1,stat="identity",fill="darkorange",color=NA) + 
    theme_classic() + theme(axis.text.x = element_blank()) + 
    geom_text_repel(size=3,min.segment.length = 0.1,segment.size = 0.2,max.overlaps = 20) + 
    geom_hline(yintercept = 0) + labs(x="TF Motifs",y="Regulation Score")
  
  gAll
}


### FINISH THIS ###
plotfigRHeatmap.2 <- function(figR.d,
                              score.cut=1,
                              DORCs=NULL,
                              TFs=NULL,
                              ... # Additional params passed to ComplexHeatmap
){
  
  
  DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
  TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()
  
  
  if(!is.null(DORCs)){
    if(!all(DORCs %in% figR.d$DORC))
      stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    DORCsToKeep <- intersect(DORCsToKeep,DORCs)
    TFsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  
  
  if(!is.null(TFs)){
    if(!all(TFs %in% figR.d$Motif))
      stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    TFsToKeep <- intersect(TFsToKeep,TFs)
    DORCsToKeep <- figR.d %>% filter(abs(Score) >= score.cut & Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
  
  
  net.d <- figR.d %>% filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>% 
    reshape2::dcast(DORC ~ Motif) %>% 
    tibble::column_to_rownames("DORC") %>% as.matrix()
  
  message("Plotting ",nrow(net.d)," DORCs x ",ncol(net.d), "TFs\n")
  
  # Heatmap view
  
  #myCols <- colorRamp2(seq(-2,2,length.out = 9),colors = jdb_palette("solar_flare"))
  myCols <- colorRampPalette(c("blue", "white", "red"))(100)
  myHeat <- Heatmap(net.d,
                    col=myCols,
                    clustering_distance_rows = "pearson",
                    clustering_distance_columns = "pearson",
                    name="Score",border = TRUE,
                    row_names_gp = gpar(fontsize=5,fontface="italic"),...)
  
  myHeat
  
}

plotfigRNetwork <- function(figR.d,
                            score.cut=1,
                            DORCs=NULL,
                            TFs=NULL,
                            size_dorc = 30,
                            size_TF = 20,
                            charge = -35,
                            weight.edges=FALSE){
  # Network view
  
  # Filter
  net.dat <-  figR.d %>% filter(abs(Score) >= score.cut)
  
  if(!is.null(DORCs))
    net.dat <- net.dat %>% filter(DORC %in% DORCs)
  
  if(!is.null(TFs))
    net.dat <- net.dat %>% filter(Motif %in% TFs)
  
  net.dat$Motif <- paste0(net.dat$Motif, ".")
  net.dat$DORC <- paste0(net.dat$DORC)
  
  
  dorcs <- data.frame(name = unique(net.dat$DORC), group = "DORC", size = size_dorc)
  tfs <- data.frame(name = unique(net.dat$Motif), group = "TF", size = size_TF)
  nodes <- rbind(dorcs,tfs)
  
  edges <- as.data.frame(net.dat)
  
  # Make edges into links (subtract 1 for 0 indexing)
  links <- data.frame(source=unlist(lapply(edges$Motif, function(x) {which(nodes$name==x)-1})), 
                      target=unlist(lapply(edges$DORC, function(x) {which(nodes$name==x)-1})), 
                      corr=edges$Corr,
                      enrichment=edges$Enrichment.P)
  
  links$Value <- scales::rescale(abs(edges$Score))*10+1
  
  # Set of colors you can choose from for TF/DORC nodes
  colors <- c("Red", "Orange", "Yellow", "Green", "Blue", "Purple", "Tomato", "Forest Green", "Sky Blue","Gray","Steelblue3","Firebrick2","Brown", "darkgrey")
  nodeColorMap <- data.frame(color = colors, hex = gplots::col2hex(colors))
  
  getColors <- function(tfColor, dorcColor = NULL) {
    temp <- c(as.character(nodeColorMap[nodeColorMap$color==tfColor,]$hex),
              as.character(nodeColorMap[nodeColorMap$color==dorcColor,]$hex))
    if (is.null(dorcColor)) {
      temp <- temp[1]
    }
    colors <- paste(temp, collapse = '", "')
    colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
    colorJS
  }
  
  forceNetwork(Links = links, 
               Nodes = nodes,
               Source = "target",
               Target = "source",
               NodeID = "name",
               #NodeID="myLabel",
               Group = "group",
               Value = "Value",
               Nodesize = "size",
               #linkWidth = 1.5, # Fixed link weight
               #linkDistance = JS("function(d) { return 5*d.value; }"),
               radiusCalculation = "Math.sqrt(d.nodesize)*2",
               arrows = FALSE,
               opacityNoHover = 0.6,
               opacity = 1,
               zoom = TRUE,
               bounded = TRUE,
               charge = charge, 
               fontSize = 13,
               legend = TRUE,
               fontFamily = "Helvetica",
               colourScale = getColors(tfColor = "Tomato",dorcColor =  "darkgrey"), # TF then DORC
               linkColour = ifelse(links$corr > 0, as.character(nodeColorMap[nodeColorMap$color=="Brown",]$hex),
                                   as.character(nodeColorMap[nodeColorMap$color=="Blue",]$hex)))
  
  
}

# Function to make J plot of significant peak-gene assocoations to call DORCs using
dorcJplot <- function(linkTab, # table returned from Links() function
                      cutoff=7, 
                      labelTop=25,
                      returnGeneList=FALSE, # Returns genes passing numPeak filter
                      cleanLabels=TRUE,
                      labelSize=4,
                      ... # Additional params passed to ggrepel
){
  
  stopifnot(all(c("peak","gene") %in% colnames(linkTab)))
  
  # Count the number of significant peak associations for each gene (without pre-filtering genes)
  numLinks <- linkTab  %>% group_by(gene) %>% tally() %>% arrange(desc(n))
  numLinks$Index <- 1:nrow(numLinks) # Add order index
  numLinks %>% as.data.frame(stringsAsFactors=FALSE) -> numLinks
  rownames(numLinks) <- numLinks$gene
  
  dorcGenes <- numLinks$gene[numLinks$n >= cutoff]
  
  numLinks <- numLinks %>%
    mutate(isDORC=ifelse(gene %in% dorcGenes,"Yes","No")) %>%
    mutate(Label=ifelse(gene %in% dorcGenes[1:labelTop],gene,""))
  
  # Plot
  dorcG <- ggplot(numLinks,aes(x=Index,y=n,color=isDORC,label=Label)) +
    geom_hline(linetype="dotted",yintercept = cutoff)+
    geom_vline(linetype="dotted",xintercept = max(numLinks[numLinks$gene %in% dorcGenes,"Index"]))+
    geom_point(size=0.8) +
    geom_line()+
    scale_color_manual(values=c("gray65","firebrick"))+
    scale_y_continuous(breaks = scales::pretty_breaks())+
    theme_classic() +
    labs(y="Number of correlated peaks",x="Ranked genes",title=paste0("Number of DORCs: ( n >= ",cutoff,") = ",length(dorcGenes)))+
    theme(axis.text = element_text(color = "black"),legend.position = "none",plot.title=element_text(hjust=0.5)) +
    scale_x_reverse() # flip so we can add labels later, if needed, with more space
  
  if(cleanLabels){
    dorcG <- dorcG + ggrepel::geom_label_repel(size=labelSize,max.iter = 1000, max.overlaps = Inf,fontface="italic",...)
  } else {
    dorcG <- dorcG + ggplot2::geom_text(size=labelSize,fontface="italic",...)
  }
  
  
  
  if(returnGeneList){
    return(dorcGenes)
    print(dorcG)
  }else(
    return(dorcG)
  )
  
  
  
}

Calculate_DORC_scores <- function(gene_list,link.df, peak_matrix){ 
  stopifnot(all(c("peak","gene") %in% colnames(link.df)))
  stopifnot(all(gene_list %in% link.df$gene))
  ds.matrix <- matrix(0, nrow = length(gene_list), ncol = ncol(peak_matrix))
  rownames(ds.matrix) <- gene_list
  colnames(ds.matrix) <- colnames(peak_matrix)
  
  ## might need to re-write this if gene_list goes larger
  for(i in 1:length(gene_list)){
    gene_name <- gene_list[[i]] 
    peaks <- links.df %>% 
      dplyr::filter(gene== gene_name) %>% 
      select(peak) %>% 
      unlist
    ds.matrix[i,] <- colSums(peak_matrix[peaks,])
  }
  ds.matrix <- Matrix::Matrix(ds.matrix, sparse = T)
  return(CreateAssayObject(counts = ds.matrix))
  
}