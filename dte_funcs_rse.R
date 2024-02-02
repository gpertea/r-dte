edgerPrep <- function(rse) {
  # this uses countsFromAbundance (as used with limma-voom)
  # the results for edgeR were nearly indistinguishable
  # with raw counts + offset and countsFromAbundnace,
  # so for code simplicity I use countsFromAbundance
  y <- DGEList(assays(rse)$scaledCounts)
  y <- y[filterByExpr(y),]
  y <- calcNormFactors(y)
  design <- model.matrix(~Condition, data=as.data.frame(colData(rse)))
  list(y, design)
}

edgerPrepFD <- function(rse) {
  # this uses filterByExpr(, design), which makes it more permissive!
  # => more features are left in
  y <- DGEList(assays(rse)$scaledCounts)
  design <- model.matrix(~Condition, data=as.data.frame(colData(rse)))
  y <- y[filterByExpr(y, design), ]
  y <- calcNormFactors(y)
  list(y, design)
}


edger <- function(rse) {
  out <- edgerPrep(rse)
  y <- out[[1]]; design <- out[[2]]
  y <- estimateDisp(y,design)
  fit <- glmFit(y,design)
  lrt <- glmLRT(fit)
  topTags(lrt, n=nrow(y), sort="none")[[1]]
}

edgerQL <- function(rse) {
  out <- edgerPrep(rse)
  y <- out[[1]]; design <- out[[2]]
  y <- estimateDisp(y,design)
  qlfit <- glmQLFit(y,design)
  qlft <- glmQLFTest(qlfit)
  topTags(qlft, n=nrow(y), sort="none")[[1]]
}

limma_voom <- function(rse) {
  scaledCounts <- tximport::makeCountsFromAbundance(assays(rse)$counts, assays(rse)$tpm,
              assays(rse)$length, countsFromAbundance = "lengthScaledTPM")
  y <- DGEList(scaledCounts)
  ## model design here:
  design <- model.matrix(~Dx, data=as.data.frame(colData(rse)))
  y <- y[filterByExpr(y, design), ]
  y <- calcNormFactors(y)
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

limmavoom <- function(rse) {
  out <- edgerPrep(rse)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

limmavoom.ftx <- function(rse) {
  out <- edgerPrepFD(rse)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

limmalogtpm <- function(rse) {
  #cts is tpm - jaffelab::expression_cutoff(cts) suggests 0.2
  tpm <- assays(rse)$tpm
  tpm <- tpm[rowMeans(tpm)>=0.24, ]
  design <- model.matrix(~Condition, data=as.data.frame(colData(rse)))
  fit = lmFit(log2(tpm + 0.5), design)
  eB <- eBayes(fit)
  topTable(eB, coef=2, n=nrow(tpm), sort="none")
}

limmatrend <- function(rse) {
  tpm <- assays(rse)$tpm
  tpm <- tpm[rowMeans(tpm)>=0.24, ]
  design <- model.matrix(~Condition, data=as.data.frame(colData(rse)))
  txExprs = log2(tpm + 0.5) # expression matrix
  arrayw = arrayWeights(txExprs, design)  # calculate weights
  fit <- lmFit(txExprs, design, weights=arrayw)
  eB <- eBayes(fit, trend = TRUE)
  topTable(eB, coef=2, n=nrow(tpm), sort="none")
}


deseq2 <- function(rse) {
  #samps <- as.data.frame(colData(rse))
  suppressMessages({
    cts <- round(assays(rse)$counts)
    mode(cts) <- 'integer'
    assays(rse)$counts <- cts
    dds <- DESeqDataSet(rse, ~Condition)
    keep <- rowSums(counts(dds) >= 10) >= ncol(rse)/2
    table(keep)
    dds <- dds[keep,]
    dds <- DESeq(dds, minReplicates=Inf, quiet=TRUE)
  })
  results(dds)
}

samseq <- function(cts, n.sub) {
  keep <- rowSums(cts >= 10) >= n.sub
  x <- round(cts[keep,])
  y <- factor(rep(1:2,each=n.sub))
  out <- capture.output({
    samfit <- SAMseq(x, y, resp.type="Two class unpaired", fdr.output=1)
  })
  sam.padj <- rep(1,nrow(x))
  names(sam.padj) <- rownames(x)
  idx <- as.numeric(samfit$siggenes.table$genes.up[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.up[,"q-value(%)"])
  idx <- as.numeric(samfit$siggenes.table$genes.lo[,"Gene Name"])
  sam.padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.lo[,"q-value(%)"])
  sam.padj
}
