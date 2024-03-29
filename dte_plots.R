rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)

cols <- c("DESeq2"="black",
          "edgeR"=rgb2(60,190,110),
          "edgeR.QL"=rgb2(0,148,80),
          "limmaVoom"=rgb2(204,121,167),
          "limmaVoomF"=rgb2(220,94,110), # filter by expression uses design => more permissive
          "limmaLogTpm"=rgb2(230,159,0),
          "limmaTrend"=rgb2(170,122,62)
          )

myicobra <- function(cd, n.sub=12, lvl="Gene", label=FALSE, zoom=FALSE) {
  if (zoom) {
    cd@padj$sleuth <- NULL
  }
  cp <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("fdrtpr","fdrtprcurve"),
                              thrs=c(.01,.05,.1))
  cobraplot <- prepare_data_for_plot(cp, colorscheme=cols)
  yrng <- if (lvl == "Gene") c(.45, 1) else c(.2, .55)
  xrng <- if (zoom) c(0, .2) else c(0,max(.21,max(fdrtpr(cp)$FDR)))
  plt <- plot_fdrtprcurve(cobraplot, plottype="points",
                          xaxisrange=xrng,
                          yaxisrange=yrng,
                          #yaxisrange=c(0,max(fdrtpr(cp)$TPR)),
                          stripsize=0,
                          title=paste0(lvl,"-level, n=",n.sub," vs ",n.sub))
  if (label)
    plt <- plt + geom_text(aes(label=method,color=method,vjust=-1))
  plt
}

fdrBreakdownGene <- function(padj.gene, alpha) {
  padj.gene <- padj.gene[apply(padj.gene, 1, function(x) !all(is.na(x))),]
  padj.gene$status <- ifelse(rownames(padj.gene) %in% union(dge.genes, dte.genes),
                             "dge",
                             ifelse(rownames(padj.gene) %in% dtu.genes,
                                    "dtu", "null"))
  ref.props <- prop.table(table(padj.gene$status[padj.gene$status != "dge"]))
  meths <- sort(setdiff(names(padj.gene), "status"))
  fdrs <- sapply(meths, function(m) {
    if (sum(padj.gene[[m]] < alpha, na.rm=TRUE) == 0) return(c(0,0))
    prop.table(table(sig=padj.gene[[m]] < alpha, padj.gene$status),1)[2,2:3]
  })
  cols <- c("goldenrod","skyblue")
  main <- paste0("Gene-level, n=", n.sub, " vs ", n.sub)
  par(mfrow=c(2,1))
  ylim <- c(0, max(.4, max(fdrs)))
  barplot(fdrs, col=cols, main=main, ylab="Total height = FDR", ylim=ylim)
  abline(h=alpha, lty=2)
  legend("top", rev(c("DTU","null")), fill=rev(cols), bg="white", title=c("Proportion by gene type"))
  barplot(cbind(reference=ref.props, t(t(fdrs)/colSums(fdrs))),
          col=cols, main=main, ylab="Proportion of FP")
}

fdrBreakdownTxp <- function(padj, alpha) {
  padj <- padj[apply(padj, 1, function(x) !all(is.na(x))),]
  iso.any <- iso.dtu | iso.dte | iso.dge
  padj$status <- rownames(padj) %in% names(iso.dtu)[iso.any]
  padj$gene <- txdf$GENEID[match(rownames(padj), txdf$TXNAME)]
  padj$type <- ifelse(padj$gene %in% dge.genes,
                      "dge",
                      ifelse(padj$gene %in% dte.genes,
                             "dte",
                             ifelse(padj$gene %in% dtu.genes,
                                    "dtu", "null")))
  ref.props <- prop.table(table(padj$type[!padj$status]))
  meths <- sort(setdiff(names(padj), c("status","gene","type")))
  fdrs <- sapply(meths, function(m) {
    #if (sum(padj.gene[[m]] < alpha, na.rm=TRUE) == 0) return(c(0,0,0))
    tot <- sum(padj[[m]] < alpha, na.rm=TRUE)
    table(FP=padj[[m]] < alpha & !padj$status, padj$type)[2,2:4]/tot
  })
  cols <- c("violet","goldenrod","skyblue")
  main <- paste0("Transcript-level, n=", n.sub, " vs ", n.sub)
  par(mfrow=c(2,1))
  ylim <- c(0, max(.4, max(fdrs)))
  barplot(fdrs, col=cols, main=main, ylab="Total height = FDR", ylim=ylim)
  abline(h=alpha, lty=2)
  legend("top", rev(c("DTE","DTU","null")), fill=rev(cols), bg="white", title="Proportion by gene type")
  barplot(cbind(reference=ref.props, t(t(fdrs)/colSums(fdrs))),
          col=cols, main=main, ylab="Proportion of FP")
}
