#' String Mismatch Handling
#'
#' Generate errors for mismatches between strings. Closest matches by
#' Levenshtein distance are calculated to offer suggestions in error message.
#'
#' @md
#' @param x A character vector. The object to search in.
#' @param y A character vector. The object to search for.
#' @param max_matches An integer. Maximum number of suggestions to print.
#' @return invisible NULL
#' @keywords internal
#' @import dplyr
#' @importFrom stringr str_match_all
#' @importFrom S4Vectors isEmpty

did_you_mean <- function(x, y, max_matches = 5) {
  if (!is.character(x) || !is.character(y)) {
    stop(strwrap(
      prefix = " ",
      initial = "", x = "Both x and y have to be of type character."
    ))
  }
  if (y %in% x) {
    stop(strwrap(
      prefix = " ", initial = "",
      x = "y found in x. Intended function usage is \"if(!y %in% x)
      did_you_mean(x,y)\""
    ))
  }

  # no match
  if ({
    adist(x, y) %>% as.vector() == nchar(x)
  } %>% all()) {
    stop(strwrap(
      prefix = " ", initial = "", x = paste(
        "No match for",
        paste("\"", y, "\"", sep = ""), "in",
        substitute(x) %>% as.character(.)
      )
    ))
  }

  # partial match
  pm <- adist(x, y) %>%
    match(., min(.)) %>%
    is.na(.) %>%
    `!`(.) %>%
    which(. == TRUE) %>%
    x[.]

  # substring match
  ssm <- str_match_all(x, y) %>%
    lapply(., FUN = S4Vectors::isEmpty) %>%
    unlist(.) %>%
    `!`(.) %>%
    x[.]
  ifelse(length(ssm) > 1,
         {
           m <- ssm
         },
         {
           m <- pm
         }
  )
  if (length(m) > max_matches) {
    error_message <- paste(
      length(m), "equally good matches for",
      paste("\"", y, "\"", sep = ""),
      "- Please specify."
    )
  } else if (length(m) > 1) {
    error_message <- paste("No match for ", paste("\"", y, "\"", sep = ""),
                           " - Did you mean one of:\n",
                           paste("\"", m, "\"", collapse = "\n", sep = ""),
                           sep = ""
    )
  } else {
    error_message <- paste(
      "No match for", paste("\"", y, "\"", sep = ""),
      "- Did you mean",
      paste("\"", m, "\"", collapse = "\n", sep = "")
    )
  }
  stop(error_message)
  return(invisible(NULL))
}
