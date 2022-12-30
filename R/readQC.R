#' Generate QC Summaries
#'
#' Summarize quality aspects and features of FASTQ files.
#'
#' @param file_paths A character vector of relative or absolute paths of fastq files.
#' @param qs A numeric vector of quantiles of PHRED scores to compute.
#' @param meta_data A dataframe storing metadata.
#' @param id_col Metadata column name where sample identifiers are stored.
#' @param extract_pattern A string specifying the sample id in the file name.
#' @param rm_pattern A string specifying the remainder to remove.
#' @return List of length 5.
#' @keywords MBT
#' @examples readQC("example.fastq", meta_data = meta.data, id_col = "sample_name", extract_pattern = "\\.\\d\\d\\d\\.", rm_pattern = "\\.")
#' @importFrom ShortRead qa
#' @importFrom reshape2 melt
#' @importFrom dplyr %>%
#' @importFrom stringr str_extract str_pad


readQC <- function(file_paths, qs = c(0.09, 0.25, 0.5, 0.75, 0.91),
                    meta_data, id_col, extract_pattern, rm_pattern){
  n <- length(file_paths)
  # init results
  res <- vector(mode = "list", length = n)
  quality.tiles <- rep(0,1000) %>% as.list() %>%
    lapply(.,`[[<-`,1, rep(0,40) %>% list()) %>% lapply(.,unlist)
  max.read.length <- 0
  # store only relevant info per file to reduce footprint
  for(i in 1:n){
    q <- ShortRead::qa(file_paths[i])
    res[[i]][["read.count"]] <- c(q[["readCounts"]]$read, file_paths[i])
    df <- q[["perCycle"]]$quality
    res[[i]][["quality.means"]] <-
      rowsum(df$Score * df$Count, df$Cycle) / rowsum(df$Count, df$Cycle) %>%
      data.frame()
    res[[i]][["quantiles"]] <-
      lapply(qs, function(q) by(df, df$Cycle, function(f)
        .get_quant(f$Score, f$Count, q), simplify = TRUE) %>%
          unlist() %>% as.vector())
    names(res[[i]][["quantiles"]]) <- as.character(qs)
    max.read.length <- max(max.read.length,
                           max(q@.srlist$perCycle$quality$Cycle))
    .phred <- q@.srlist$perCycle$quality[,c(1,3,4)] %>% data.frame()
    for(j in 1:nrow(.phred)){
      .cycle <- .phred[j, 1]
      .score <- .phred[j, 2]
      .count <- .phred[j, 3]
      quality.tiles[[.cycle]][[.score]] <-
        quality.tiles[[.cycle]][[.score]] + .count
    }
    # print progress
    cat('\014\n')
    extra <- nchar('||100%')
    width <- options()$width
    step <- round(i / n * (width - extra))
    text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                    strrep(' ', width - step - extra), round(i / n * 100))
    cat(text)
  }
  cat('\n')
  # aggregate qc results
  read.counts <- lapply(res, `[[`, 1) %>%
    data.frame() %>% t() %>% data.frame() %>%
    mutate(X1 = as.numeric(X1)) %>%
    `names<-`(c("Read.Counts","SampleID"))
  # add metadata
  read.counts$Sample <- stringr::str_extract(read.counts$SampleID, extract_pattern) %>%
    gsub(rm_pattern, "", .)
  read.counts <- merge(read.counts, meta_data, by.x = "Sample", by.y = id_col)
  read.means <- lapply(res, `[[`, 2) %>% lapply(., `[[`, 1) %>%
    data.frame() %>% t() %>% melt()
  read.quants <- lapply(res, `[[`, 3) %>%
    lapply(.,as.data.frame, check.names=F) %>%
    mapply(function(x) "[<-"(x, "Cycle", value = 1:nrow(x)) ,
           ., SIMPLIFY = FALSE) %>%
    mapply(melt, ., id.vars = "Cycle", SIMPLIFY = FALSE) %>%
    mapply(cbind, ., "SampleID" = 1:n, SIMPLIFY = FALSE) %>% bind_rows()
  quality.tiles <- quality.tiles[1:max.read.length]
  quality.tiles <- quality.tiles %>% as.data.frame() %>% t() %>%
    `rownames<-`(1:max.read.length) %>% melt() %>%
    `names<-`(c("Cycle", "Score", "Count"))
  quality.tiles$Count[quality.tiles$Count==0] <- NA
  # Oedenwinkel specific
  #read.quants$Sample <- stringr::str_pad(read.quants$Sample, 3, pad = "0")
  #read.quants <- merge(read.quants, meta_data, by.x = "Sample", by.y = id_col)
  res <- list(read.counts, read.means, read.quants,
               quality.tiles, max.read.length)
  names(res) <- c("read.counts", "read.means", "read.quants",
                   "quality.tiles", "max.read.length")
  return(res)
}


.get_quant <- function(score, count, q) {
  score[which(cumsum(count)/sum(count) >= q)][[1]]
}
