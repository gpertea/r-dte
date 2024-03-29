edgerPrep <- function(cts, n.sub) {
  # this uses countsFromAbundance (as used with limma-voom)
  # the results for edgeR were nearly indistinguishable
  # with raw counts + offset and countsFromAbundnace,
  # so for code simplicity I use countsFromAbundance
  y <- DGEList(cts)
  y <- y[filterByExpr(y),]
  y <- calcNormFactors(y)
  condition <- factor(rep(1:2,each=n.sub))
  design <- model.matrix(~condition)
  list(y, design)
}

edgerPrepFD <- function(cts, n.sub) {
  # this uses filterByExpr(, design), which makes it more permissive!
  # => more features are left in
  y <- DGEList(cts)
  condition <- factor(rep(1:2,each=n.sub))
  design <- model.matrix(~condition)
  y <- y[filterByExpr(y, design), ]
  y <- calcNormFactors(y)
  list(y, design)
}


edger <- function(cts, n.sub) {
  out <- edgerPrep(cts, n.sub)
  y <- out[[1]]; design <- out[[2]]
  y <- estimateDisp(y,design)
  fit <- glmFit(y,design)
  lrt <- glmLRT(fit)
  topTags(lrt, n=nrow(y), sort="none")[[1]]
}

edgerQL <- function(cts, n.sub) {
  out <- edgerPrep(cts, n.sub)
  y <- out[[1]]; design <- out[[2]]
  y <- estimateDisp(y,design)
  qlfit <- glmQLFit(y,design)
  qlft <- glmQLFTest(qlfit)
  topTags(qlft, n=nrow(y), sort="none")[[1]]
}

limmavoom <- function(cts, n.sub) {
  out <- edgerPrep(cts, n.sub)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

limmavoom.ftx <- function(cts, n.sub) {
  out <- edgerPrepFD(cts, n.sub)
  y <- out[[1]]; design <- out[[2]]
  v <- voom(y, design)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, coef=2, n=nrow(y), sort="none")
}

limmalogtpm <- function(cts, n.sub) {
  #cts is tpm - jaffelab::expression_cutoff(cts) suggests 0.43
  tpm <- cts[rowMeans(cts)>=0.43, ]
  condition <- factor(rep(1:2,each=n.sub))
  design <- model.matrix(~condition)
  fit = lmFit(log2(tpm + 0.5), design)
  eB <- eBayes(fit)
  topTable(eB, coef=2, n=nrow(tpm), sort="none")
}


limmatrend <- function(cts, n.sub) {
  #cts is tpm - jaffelab::expression_cutoff(cts) suggests 0.43
  tpm <- cts[rowMeans(cts)>=0.43, ]
  condition <- factor(rep(1:2,each=n.sub))
  design <- model.matrix(~condition)
  txExprs = log2(tpm+0.5) # expression matrix
  arrayw = arrayWeights(txExprs, design)  # calculate weights
  fit <- lmFit(txExprs, design, weights=arrayw)
  eB <- eBayes(fit, trend = TRUE)
  topTable(eB, coef=2, n=nrow(tpm), sort="none")
}


deseq2 <- function(txi, n.sub) {
  n <- ncol(txi$counts)/2
  idx <- c(1:n.sub, n + 1:n.sub)
  samps <- data.frame(condition=factor(rep(1:2,each=n)))
  suppressMessages({
    dds <- DESeqDataSetFromTximport(txi, samps, ~condition)
    dds <- dds[,idx]
    keep <- rowSums(counts(dds) >= 10) >= n.sub
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
