#####################################################

options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

suppressWarnings(suppressMessages(library("optparse")))
library("stats")

option_list = list(
  make_option(c("-c", "--countsfile"), type="character", default=NULL, 
              help="count dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.csv", 
              help="output CSV filename", metavar="character"),
  make_option(c("-p", "--out2"), type="character", default="out.pdf", 
              help="output PDF filename", metavar="character"),
  make_option(c("-w", "--wilcoxon"), action="store_true", default=FALSE, 
              help="compute Mann-Whitney-Wilcoxon test", metavar="character"),
  make_option(c("-d", "--deseqpar"),action="store_true",default=FALSE, 
              help="Run DESeq with parametric ", metavar="character"),
  make_option(c("-l", "--deseqloc"),action="store_true",default=FALSE, 
              help="Run DESeq with local", metavar="character"),
  make_option(c("-x", "--deseq2par"),action="store_true",default=FALSE, 
              help="Run DESeq with parametric ", metavar="character"),
  make_option(c("-y", "--deseq2loc"),action="store_true",default=FALSE, 
              help="Run DESeq with local", metavar="character"),
  make_option(c("-m", "--metagenomeseqLog"), action="store_true",default=FALSE, 
              help="run metagenomeSeq with LogNormal model", metavar="character"),
  make_option(c("-z", "--metagenomeseqZIG"), action="store_true",default=FALSE, 
              help="run metagenomeSeq with ZIG model", metavar="character"),
  make_option(c("-e", "--edgeR"), action="store_true",default=FALSE, 
              help="run edgeR ", metavar="character"),
  make_option(c("-f", "--filter"), type="numeric", default=0,
              help="filter out genes with less counts in one of conditions", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=TRUE);

datafile <- opt$options$countsfile
print(paste("Running differential expression on file:", datafile))

args <- opt$args
countTable = read.table( datafile, sep="\t", header=TRUE, row.names=1 )

if ((length(args) != length(colnames(countTable))) || (sum(args == "0") + sum(args == "1") + sum(args == '-') != length(args))) {
  cat(" You need to specify if sample is treated (1) or untreated(0) for every sample\n")
  cat(args)
  q(save="no", status=1)
}


############################### Prepare data ##############################################

countTable[,] <- round(countTable[,])
countTable[, args == '-'] <- NULL

f <- function(x) {
  f <- ifelse(x=="0", "untreated", "treated")
}
condition <- sapply(args[args != '-'], f)
n <- length(condition)

print("Treated samples:")
print(colnames(countTable[condition == "treated"]))
print("Untreated samples:")
print(colnames(countTable[condition == "untreated"]))

n1<- sum(condition == 'untreated')
n2<-sum(condition == 'treated')

sums1 <- rowSums(countTable[, condition == "untreated", drop=FALSE ])
sums2 <- rowSums(countTable[, condition == "treated", drop=FALSE ])

data <- countTable
data$name <- rownames(data)
data$mean <- (sums1+sums2)/n
data$mean1 <- sums1/n1
data$mean2 <- sums2/n2
data$fold <- data$mean1/data$mean2
data$log2fold <- sapply(data$fold, log2)

NON_PV_FIELDS=6

MINREADS <- 1
NMINREADS <- n-1

successfull <- 0    #how many tests run successfully

########### remove some reads with not enough reads in condition ####################
if (opt$options$filter != 0){
  MINREADSCOND <- opt$options$filter
  to_remove <- (sums2/n2<=MINREADSCOND) | (sums1/n1<=MINREADSCOND) 
  print(paste("Removing", sum(to_remove), "genes with <=", MINREADSCOND, "counts in condition"))
  
  #int <- data[to_remove,]
  #write.table(int[order(abs(int$log2fold)),(n+1):dim(int)[2]], paste(opt$options$out, '/to_remove.txt', sep=''), quote=FALSE, sep='\t',
  #                        col.names = FALSE)
  countTable <- countTable[!to_remove,] 
  #data <- data[!to_remove,]
}

########## run simpler test first - Mann-Whitney-Wilcoxon, H0: equal means, no assumption on distr

run.wilcoxon <- function(s, i) {
  pv <- NA
  pv <- tryCatch({
      wilcox.test(as.numeric(s[i,condition=="untreated"]), 
              as.numeric(s[i,condition=="treated"]), 
              alternative = "two.sided")$p.value
    }, warning = function(w) {
      NA
  })
  return(pv)
}


if (opt$options$wilcoxon) {
  result = tryCatch({
    print("Running Mann-Whitney-Wilcoxon test")
    pv = rep(NA,dim(countTable)[1])
    names(pv) <- rownames(countTable)
    normfac <- colSums(countTable)
    s <- scale(countTable, center=FALSE, scale=normfac/1000000)
    
    for (i in 1:dim(countTable)[1]) { pv[i]=run.wilcoxon(s, i)  }    
    
    data$MWpval <- sapply(rownames(data), function(x) pv[x])
    
    #print(head(pv))
    #data$MWpadj <- p.adjust(pv, method="BH")
    rm(s)    
    successfull<-successfull+1
    print("Mann-Whitney-Wilcoxon test finished successfully.")
  }, warning = function(w) {
    print(w)
  }, error = function(e) {
    print("Error occured running Wilcoxon test")
    print(e)
  })
}

####################### edgeR #############################
if (opt$options$edgeR) {
  result = tryCatch({
    print("Running edgeR")
    suppressWarnings(suppressMessages(library(edgeR)))
    group <- sapply(condition, function(x) ifelse(x=="treated", 1,2))
    y <- DGEList(counts=countTable,group=group)
    y <- calcNormFactors(y)
    keep <- rowSums(cpm(y)>=MINREADS) >= NMINREADS      #filter?
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    
    design <- model.matrix(~group)
    y <- estimateDisp(y,design)
    
    #To perform quasi-likelihood F-tests:   for small number of replicates
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    topTags(qlf)
    
    data$edgeRpv <- sapply(rownames(data), function(x) qlf$table[x,"PValue"])
    #padj <- p.adjust(qlf$table[,"PValue"]) 
    #data$edgeRpadj <- sapply(rownames(data), function(x) padj[x])
    
    #To perform likelihood ratio tests:
    #fit <- glmFit(y,design)
    #lrt <- glmLRT(fit,coef=2)
    #topTags(lrt)
    successfull<-successfull+1
    print("edgeR finished successfully.")
  }, warning = function(w) {
    print("edgeR warning:")
    print(w)
  }, error = function(e) {
    print("Error occured running edgeR")
    print(e)
  })
  rm(y, fit, qlf)
}

########################### DESeq or DESeq2 plot output #####################
if (opt$options$deseqpar | opt$options$deseqloc | opt$options$deseq2par | opt$options$deseq2loc) {
  title=paste("DESeq results"  ) 
  pdf(opt$options$out2, paper='a4', onefile=TRUE, title=title)
}

deseq.generalinfo <- function(cds, title) {
  suppressWarnings(suppressMessages(library(gplots)))
  par(mar=c(0,0,0,0))
  tab <- data.frame(names=colnames(countTable))
  #for my data with long colnames:
  #y <- strsplit(colnames(countTable),".",fixed=TRUE)
  #tab$names <- lapply(y, function(x) {paste(x[6],x[10],sep=".")})
  
  tab$names <- colnames(countTable)
  tab$condition <- condition
  tab$sizeFactors <-  lapply(sizeFactors(cds), function(x) { formatC(x, digits=4, format="g", flag="-") })
  tab$nReads <- colSums(countTable)
  colnames(tab) <- c("sample", "condition", "size factors", "no. reads")
  
  textplot(tab, cex=1.5, show.rownames=F)
  mtext(title, outer=TRUE, cex=1.3, line=-1) 
}

plot.deseq1 <- function(cds, fitType, title) {
  par(oma=c(0,0,3,0), mar=c(5,5,3,5))
  #fitInfo( cds ) 
  DESeq::plotDispEsts( cds, main="Empirical (black dots) and fitted (red lines) dispersion values \nplotted against the mean of the normalised counts. " )
    
  mtext(title, outer=TRUE, cex=1.5, line=1) 
  res =  DESeq::nbinomTest( cds, "untreated", "treated" )
  DESeq::plotMA(res, main="Plot of normalised mean versus log2 fold change \nfor the contrast untreated versus treated")    
    
  hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="Histogram of p-values")
  par(xaxs='i', yaxs='i' )
    
  #res <- res[ order(res$pval), ]
  #write.table(res, paste(outdir, '/', fitType, '.csv', sep=''), append=FALSE, quote=FALSE, row.names = FALSE, sep="\t")
  return(res)
}

plot.deseq2 <- function(cds, fitType, title) {
  par(oma=c(0,0,3,0), mar=c(5,5,3,5))
  DESeq2::plotDispEsts( cds, main="Empirical (black dots) and fitted (red lines) dispersion values \nplotted against the mean of the normalised counts. " )
  
  mtext(title, outer=TRUE, cex=1.5, line=1) 
  #res = DESeq2::nbinomLRT( cds )
  cds = DESeq2::nbinomWaldTest( cds )    
  res <- results(cds)
  DESeq2::plotMA(res, main="Plot of normalised mean versus log2 fold change \nfor the contrast untreated versus treated")

  hist(res$pvalue, breaks=100, col="skyblue", border="slateblue", main="Histogram of p-values")
  par(xaxs='i', yaxs='i' )
    
  #res <- res[ order(res$pvalue), ] 
  #write.table(res, paste(outdir, '/', fitType, '.csv', sep=''), append=FALSE, quote=FALSE, row.names = FALSE, sep="\t")
  return(res)
}

#######  Run DESeq ###########
if (opt$options$deseqpar | opt$options$deseqloc) {
  result = tryCatch({
    suppressWarnings(suppressMessages(library( "DESeq" )))

    cds = DESeq::newCountDataSet( countTable, condition )
    cds = DESeq::estimateSizeFactors( cds )
    print(paste("Running DESeq"))
    
    if ( opt$options$deseqpar ) {
      tryCatch({
        print(paste("Running DESeq parametric"))
        title=paste("DESeq on", nrow(countTable) , "genes.", "Fit type: parametric")
        deseq.generalinfo(cds, title)
        
        cdsp <- DESeq::estimateDispersions( cds, fitType="parametric" )
        res <- plot.deseq1(cdsp, "parametric", title)
        data$DESeqparpv <- sapply(rownames(data), function(x) ifelse(sum(res$id==x) == 1, res[res$id==x,]$pval, NA))
        #data$DESeqparpadj <- sapply(rownames(data), function(x) res[res$id==x,]$padj)
        successfull<-successfull+1
        print("DESeq parametric finished successfully.")
      }, warning = function(w) {
        print("DESeq parametric warning:"); print(w)
     }, error = function(e) {
        print("Error occured running DESeq parametric -> running local instead"); print(e)
        opt$options$deseqloc <- TRUE
      }    )
    } 
    
    if ( opt$options$deseqloc ) {
      tryCatch({
        print(paste("Running DESeq local"))
        title=paste("DESeq on", nrow(countTable) , "genes.", "Fit type: local")
        deseq.generalinfo(cds, title)
      
        cdsl <- DESeq::estimateDispersions( cds, method="blind", sharingMode="fit-only", fitType="local" )
        res <- plot.deseq1(cdsl, "local", title)
        data$DESeqlocpv <- sapply(rownames(data), function(x) ifelse(sum(res$id==x) == 1, res[res$id==x,]$pval, NA))
        #data$DESeqlocpadj <- sapply(rownames(data), function(x) res[res$id==x,]$padj)
        successfull<-successfull+1
        print("DESeq local finished successfully.")    
      }, warning = function(w) {
        print("DESeq local warning:"); print(w)
     }, error = function(e) {
        print("Error occured running DESeq local"); print(e)
        opt$options$deseqloc <- TRUE
     })
    }
  }, warning = function(w) { print("DESeq warning:"); print(w)
  }, error = function(e) { print("Error occured running DESeq"); print(e)
  })
  rm(cds)
  garbage <- dev.off()
}

################### DESeq2 ##############################

if (opt$options$deseq2par | opt$options$deseq2loc) {
  result = tryCatch({
    suppressWarnings(suppressMessages(library( "DESeq2" )))
    print(paste("Running DESeq2"))

    coldata = data.frame(colnames(countTable), condition)
    countData <- data.frame(countTable)
    colnames(countData) <- NULL
    colnames(coldata)=c("sample", "condition")
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = coldata, design = ~ condition)    
    dds = DESeq2::estimateSizeFactors( dds )
    
    if ( opt$options$deseq2par ) {
      tryCatch({
        print(paste("Running DESeq2 parametric"))
        title=paste("DESeq2 on", nrow(countTable) , "genes.", "Fit type: parametric")
        deseq.generalinfo(dds, title)
        
        ddsp <- DESeq2::estimateDispersions( dds, fitType="parametric" )
        res <- plot.deseq2(ddsp, "parametric", title)
        data$DESeq2parpv <- sapply(rownames(data), function(x) ifelse(x %in% rownames(res), res[x,]$pvalue, NA))
        
        #data$DESeq2parpadj <- sapply(rownames(data), function(x) res[res$id==x,]$padj)
        successfull<-successfull+1
        print("DESeq2 parametric finished successfully.")
      }, error = function(e) {
        print("Error occured running DESeq2 parametric -> running local instead")
        print(e)
        opt$options$deseq2loc <- TRUE
      }    )
    } 
    
    if ( opt$options$deseq2loc ) {
      print(paste("Running DESeq2 local"))        
      title=paste("DESeq2 on", nrow(countTable) , "genes.", "Fit type: local")
      deseq.generalinfo(dds, title)
      
      ddsl <- DESeq2::estimateDispersions( dds, fitType="local" )
      res <- plot.deseq2(ddsl, "local", title )
      data$DESeq2locpv <- sapply(rownames(data), function(x) ifelse(x %in% rownames(res), res[x,]$pvalue, NA))
      
      #data$DESeq2locpadj <- sapply(rownames(data), function(x) res[res$id==x,]$padj)
      successfull<-successfull+1
      print("DESeq2 local finished successfully.")
    }  
  }, warning = function(w) {
    print("DESeq2 warning:")
    print(w)
  }, error = function(e) {
    print("Error occured running DESeq2")
    print(e)
  } )
  rm(countData, dds)
  garbage <- dev.off()
}

################## metagenomeSeq #######################
if ( opt$options$metagenomeseqLog | opt$options$metagenomeseqZIG ) {
  result = tryCatch({
    suppressWarnings(suppressMessages(library(metagenomeSeq)))
    print("Running MetagenomeSeq")
    metaData <- loadMeta(datafile, sep = '\t' )
    #metaData <- load_meta(datafile,sep = '\t' )
    mrexpr <- newMRexperiment(metaData$counts)
    samplesToKeep = which(args != '-')
    mrexpr = mrexpr[, samplesToKeep]
    #filtrowanie
    print(paste("MetagenomeSeq: genes before filtering", dim(MRcounts(mrexpr))[1]))
    
    p = cumNormStatFast(mrexpr)
    mrexpr = cumNorm(mrexpr, p = p)
    mrexpr = filterData(mrexpr, present = NMINREADS, depth = MINREADS)
    print(paste("MetagenomeSeq: genes after filtering", dim(MRcounts(mrexpr))[1]))
    
    #### another filtering?
    #rareFeatures = which(rowSums(MRcounts(mrexpr, norm=T) < 1) > 1)
    #mrexpr = mrexpr[-rareFeatures, ]
    #print(paste("MetagenomeSeq filtered out", length(rareFeatures), "reads leaving", dim(MRcounts(mrexpr))[1]))
    
    ###export counts
    #mat = MRcounts(mrexpr, norm = TRUE, log = TRUE)
    #exportMat(mat, file = outputfile6)
    #exportStats(mrexpr, file= outputfile6)
    
    if ( opt$options$metagenomeseqLog ) {      ## log-normal model
      result = tryCatch({
        pd <- pData(mrexpr)
        mod <- model.matrix(~condition, data = pd)
        fitModel = fitFeatureModel(mrexpr, mod)
        #fitModel = fitLogNormal(mrexpr, mod)
        head(MRcoefs(fitModel)$pvalues)
        
        data$metagenomeSeqLogNorm.pv <- sapply(rownames(data), function(x) fitModel$pvalues[x])
        #padj <- p.adjust(fitMode$pvalues, method="BH")
        #data$metagenomeSeqLogNorm.padj <- sapply(rownames(data), function(x) padj[x])        
        successfull<-successfull+1
        print("MetagenomeSeq LogNorm finished successfully.")
      }, error = function(e) {
        print("Error occured running MetagenomeSeq LogNorm")
        print(e)
      } )
    }
    if ( opt$options$metagenomeseqZIG ) {      ##Zero-Inflated-Gaussian model 
      result = tryCatch({
        normFactor = normFactors(mrexpr)
        normFactor = log2(normFactor/median(normFactor) + 1)
        mod2 = model.matrix(~condition + normFactor)
        settings = zigControl(maxit = 10, verbose = TRUE)
        fitZigModel = fitZig(obj = mrexpr, mod = mod2, useCSSoffset = FALSE,
                             control = settings)
        
        data$metagenomeSeqZIG.pv <- sapply(rownames(data), function(x) fitZigModel$eb$p.value[,2][x])
        #padj <- p.adjust(fitMode$pvalues, method="BH")
        #data$metagenomeSeqZIG.padj <- sapply(rownames(data), function(x) padj[x])
        successfull<-successfull+1
        print("MetagenomeSeq ZIG finished successfully.")
      }, error = function(e) {
        print("Error occured running MetagenomeSeq ZIG")
        print(e)
      } )
    }  
  }, error = function(e) {
    print("Error occured preprocessing MetagenomeSeq")
    print(e)
  } )
}

################# WRITE SUMMARY
#print(warnings())

if (successfull == 0) {
  write("No method finished successfully", stderr())
  q(save="no", status=1)
}

data <- data[, (n+1):dim(data)[2]]
last_pv <- dim(data)[2]

if (successfull > 1) {
  data$toCheck <- apply(data[,(NON_PV_FIELDS+1):last_pv], 1, function(x) sum(x < 0.05, na.rm = T))
  print(paste("Number of genes with p-value < 0.05 according to at least 2 methods:", sum(data$toCheck >= 2, na.rm = TRUE)))
  write.table(data[order(data$toCheck,  decreasing = T),], opt$options$out, quote=FALSE, sep='\t', row.names = F,
              col.names = TRUE)
} else {
  print(paste("Number of genes with p-value < 0.05:", sum(data[,last_pv] < 0.05, na.rm = TRUE)))  
  write.table(data[order(data[,last_pv],  decreasing = F),], opt$options$out, quote=FALSE, sep='\t', row.names = F,
              col.names = TRUE)
}

q(save="no", status=0)
