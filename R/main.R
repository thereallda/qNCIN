#' Generate scale factors based on spike-ins
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for fly spike-in id, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#'
#' @return Vector of scale factors
#' @export
#'
#' @importFrom edgeR calcNormFactors
calcScaleFactor <- function(data,
                            spike.in.prefix = NULL,
                            enrich.group,
                            input.id = "Input",
                            enrich.id = "Enrich") {
  #TODO: store parameters in object
  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from spike-in
  counts_spk <- data[grep(spike.in.prefix, rownames(data)),]

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # compute scale factors using TMM method
  nf1 <- edgeR::calcNormFactors(counts_spk, refColumn = input.idx)
  sf1 <- nf1/colSums(counts_spk)*1e6

  # adjust scale factor by ratio of spike-in/non-spike-in
  lib_ratio <- colSums(counts_spk)/colSums(data)
  lib_ratio <- lib_ratio/(1-lib_ratio)
  lib_ratio <- c(lib_ratio[input.idx]/mean(lib_ratio[input.idx]),
                 lib_ratio[enrich.idx]/mean(lib_ratio[enrich.idx]))
  scale.factor <- sf1*lib_ratio

  return(scale.factor)
}

#' Generate adjust factors based on non-specific enriched genes
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for fly spike-in id, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#' @param scale.factor Vector of scale factors generated from \code{calcScaleFactor}.
#' @param prop.top.enrich Proportion of top-enriched genes to used for
#' adjustment of non-specific enrichment, default used all genes.
#' @param top.down Whether using top-to-down or down-to-top enriched gens for
#' the adjustment of non-specific enrichment, default: TRUE (top-down)
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return Vector of adjust factors
#' @export
#'
#' @importFrom stats lm
#' @importFrom utils head
calcAdjustFactor <- function(data,
                             spike.in.prefix = NULL,
                             enrich.group,
                             input.id = "Input",
                             enrich.id = "Enrich",
                             scale.factor,
                             prop.top.enrich = 1,
                             top.down = TRUE,
                             pseudo.count = 1) {
  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # # ensure each samples have at least one count
  # counts_sam <- counts_sam[rowSums(counts_sam > 0) == ncol(counts_sam),]

  # add pseudo counts to shrink ratio
  counts_sam <- addPseudoCount(counts_sam, pseudo.count = pseudo.count)

  # scale sample counts by scale factor
  counts_sam_scale <- t(t(counts_sam)*scale.factor)

  # calculate empirical NCIN ratio
  eratio_df <- counts_sam_scale[,enrich.idx]/counts_sam_scale[,input.idx]

  # get all non-specific enriched genes based on empirical ratio
  neg_all <- eratio_df[rowSums(eratio_df>1) > 0.7*ncol(eratio_df),]
  neg_avg <- rowMeans(neg_all)
  neg_all <- cbind(neg_all, neg_avg)

  # order neg based on their average in decreasing
  neg_all <- neg_all[order(neg_all[,"neg_avg"], decreasing = top.down),]

  # used prop.top.enrich of non-specific enriched genes for adjustment
  neg_all <- utils::head(neg_all, n = floor(nrow(neg_all)*prop.top.enrich))

  # create centroid based on non-specific enriched genes from input samples
  centroid <- rowMeans(counts_sam_scale[rownames(neg_all), input.idx])

  # linear regression between centroid and each sample
  reg_df <- log2(as.data.frame(cbind(counts_sam_scale[rownames(neg_all),], centroid)))
  adjust.factor <- c()
  for (id in colnames(counts_sam)) {
    dfi <- reg_df[,c("centroid", id)]
    colnames(dfi) <- c("centroid", "x")
    lmi <- stats::lm(centroid ~ 1 + offset(x), data = dfi)
    adjust.factor[id] <- 2^lmi$coefficients
  }

  return(adjust.factor)
}

#' Compute NCIN ratio
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the
#'   number of samples and p is the number of features.
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g.,
#'  "^FB" stands for fly spike-in id, default: NULL.
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input".
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".
#' @param scale.factor Vector of scale factors generated from \code{calcScaleFactor}.
#' @param adjust.factor Vector of adjust factors generated from \code{calcAdjustFactor}.
#' @param filter Whether to keep genes with ratio smaller than one in all samples, default: TRUE
#' @param pseudo.count A numeric scalar of pseudo-counts to be added to each gene, default: 1.
#'
#' @return data.frame containing the ratio of each genes from each samples,
#' and the average and standard deviation of computed ratio from each gene.
#' @export
#'
#' @importFrom stats sd
calcNCIN <- function(data,
                     spike.in.prefix = NULL,
                     enrich.group,
                     input.id = "Input",
                     enrich.id = "Enrich",
                     scale.factor,
                     adjust.factor,
                     filter = TRUE,
                     pseudo.count = 1) {

  # initialize parameters
  input.idx <- grep(input.id, enrich.group)
  enrich.idx <- grep(enrich.id, enrich.group)

  # get counts from sample
  counts_sam <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]

  # # ensure each samples have at least one count
  # counts_sam <- counts_sam[rowSums(counts_sam > 0) == ncol(counts_sam),]

  # add pseudo counts to shrink ratio
  counts_sam <- addPseudoCount(counts_sam, pseudo.count = pseudo.count)

  # scale sample counts by scale factor
  counts_sam_scale <- t(t(counts_sam)*scale.factor)

  # further adjust sample counts by adjust factor
  counts_sam_adj <- t(t(counts_sam_scale)*adjust.factor)

  # calculate quantitative ratio
  ratio.df <- counts_sam_adj[,enrich.idx]/counts_sam_adj[,input.idx]

  # keep gene with ratio smaller than 1 in all samples
  if (filter) {
    ratio.df <- ratio.df[rowSums(ratio.df<1) == ncol(ratio.df), ]
  }

  # return summarized statistics
  ratio.avg <- rowMeans(ratio.df)
  ratio.sd <- apply(ratio.df, 1, stats::sd)
  ratio.df <- as.data.frame(cbind(ratio.df, ratio.avg, ratio.sd))

  return(ratio.df)
}
