#' Distance Metric
#'
#' Compute distances of the feature table of a _phyloseq_ object.
#'
#' @param obj A phyloseq object.
#' @param dis Name of distance to use. One of _philr_, _u.unifrac_, _w.unifrac_,
#' or a metric accepted by the _method_ argument in *vegan::vegdist()*.
#' @param part.weights See *philr::philr()*.
#' @param ilr.weights See *philr::philr()*.
#' @param pseudocount Passed to *vegan::vegdist()* if dis is _aitchison_ or
#' _robust.aitchison_, or added to ALL counts of feauture table if
#' _dis = philr_.
#' @param normalized See *phyloseq::UniFrac()*.
#' @param parallel See *phyloseq::UniFrac()*.
#' @param ... Passed to *vegan::vegdist()* and *philr::philr()*.
#' @return distance matrix
#' @keywords MBT
#' @export
#' @examples
#' mbt_dist <- MBT::dist(GlobalPatterns, dis = w.unifrac)
#' @import ape vegan dplyr phyloseq philr

dist <- function(obj, dis = "bray", part.weights = "uniform",
                    ilr.weights = "uniform", pseudocount = 1e-12,
                    normalized = TRUE, parallel = FALSE, ...)
{
  tar <- obj@otu_table@taxa_are_rows
  ft <- otu_table(obj)
  if(tar){
    ft <- t(ft)
  }
  ft_bool <- ifelse(ft > 0, 1, 0)
  x <- rowSums(ft_bool) %>% data.frame()
  if(any(x==0)) {
    warning("One or more samples contain zero features and will be removed.")
  }
  ft <- ft[x > 1, ] #should be 0?
  if (dis %in% c("aitchison", "robust.aitchison")) {
    dist <- vegan::vegdist(ft, method = dis, pseudocount = pseudocount, ...)
  } else if(!dis %in% c("u.unifrac", "w.unifrac", "philr")) {
    dist <- vegan::vegdist(ft, method = dis, ...)
  } else {
    if(dis != "philr") {
      dist <- phyloseq::UniFrac(obj, weighted = (dis == "w.unifrac"),
                                normalized = normalized, parallel = parallel)
    } else {
      obj@otu_table <- otu_table(obj) + pseudocount
      philr <- philr::philr(obj, part.weights = part.weights,
                            ilr.weights = ilr.weights, ...)
      dist <- vegan::vegdist(philr, method = "euclidian")
    }
  }
  return(dist)
}
