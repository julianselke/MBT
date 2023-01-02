#' Show Loop Progress
#'
#' Print progress bar, percentage, and current iteration to console.
#'
#' @param i An integer. The loop counter.
#' @param n An integer. The maximum loop count.
#' @param names A character vector. The names of the dimension being looped
#' over, e.g. row.names.
#' @keywords internal

show_progress <- function(i, n, names) {
  cat('\014\n')
  extra <- nchar('||100%')
  width <- options()$width
  step <- round(i / n * (width - extra))
  text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                  strrep(' ', width - step - extra), round(i / n * 100))
  cat(text)
  cat("\nProcessing ", names[i])
}
