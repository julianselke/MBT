#' Calculate Distance/Dissimilarity Metric
#'
#' Compute a distance matrix from the feature table of a phyloseq object.
#'
#' **SAMPLES ARE ROWS**
#'
#' @md
#' @param obj A phyloseq object.
#' @param dis Name of distance to use. One of philr, u.unifrac, w.unifrac,
#' or a metric accepted by the method argument in \link[vegan]{vegdist}.
#' @param part.weights See \link[philr]{philr}.
#' @param ilr.weights See \link[philr]{philr}.
#' @param pseudocount Passed to \link[vegan]{vegdist} if dis is aitchison or
#' robust.aitchison, or added to ALL counts of the feature table if
#' dis = philr.
#' @param normalized See \link[phyloseq]{UniFrac}.
#' @param parallel See \link[phyloseq]{UniFrac}.
#' @param ... Passed to \link[vegan]{vegdist} and \link[philr]{philr}.
#' @return A distance matrix.
#' @keywords internal
#' @import ape
#' @importFrom vegan vegdist
#' @import dplyr
#' @import phyloseq
#' @importFrom philr philr

dist <- function(obj,
                 dis = "bray",
                 part.weights = "uniform",
                 ilr.weights = "uniform",
                 pseudocount = 1e-12,
                 normalized = TRUE,
                 parallel = FALSE,
                 ...) {
  # Check orientation of feature table in phyloseq object since there is
  # no internal control for assumption "samples are rows" in philr.
  # philr::philr --> philr.phyloseq --> philr.data.frame
  if (isTRUE(obj@otu_table@taxa_are_rows)) obj <- phyloseq::t(obj)
  ft <- otu_table(obj)
  ft_bool <- ifelse(ft > 0, 1, 0)
  x <- rowSums(ft_bool) %>% data.frame()
  if (any(x == 0)) {
    warning("One or more samples contain zero features and will be removed.")
  }
  ft <- ft[x > 0, ]
  if (dis %in% c("aitchison", "robust.aitchison")) {
    dist <- vegan::vegdist(ft, method = dis, pseudocount = pseudocount, ...)
  } else if (!dis %in% c("u.unifrac", "w.unifrac", "philr")) {
    dist <- vegan::vegdist(ft, method = dis, ...)
  } else {
    if (dis != "philr") {
      dist <- phyloseq::UniFrac(obj,
                                weighted = (dis == "w.unifrac"),
                                normalized = normalized,
                                parallel = parallel
      )
    } else {
      otu_table(obj) <- otu_table(obj) + pseudocount
      ph.ft <- philr::philr(otu_table(obj),
                            phy_tree(obj),
                            part.weights = part.weights,
                            ilr.weights = ilr.weights,
                            ...
      )
      dist <- vegan::vegdist(ph.ft, method = "euclidian")
    }
  }
  return(dist)
}
