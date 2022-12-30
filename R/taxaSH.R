#' Taxa Overview
#'
#' Start a *Shiny App* to explore taxonomic composition of a _phyloseq_ dataset.
#'
#' @param obj A _phyloseq_ object.
#' @param plot_height An integer number specifying the plot height in *px*.
#' Defaults to 1.6 times the size of a default plot window opened from Rscript.
#' @keywords MBT
#' @export
#' @examples
#' taxaSH(GlobalPatterns)
#' @import shinyWidgets plotly purrr shiny dplyr ggplot2 FacileViz phyloseq reshape2


taxaSH <- function (obj, plot_height = NA){
  if(!is.na(plot_height)){
    if(is.numeric(plot_height) == FALSE | plot_height <= 0){
      stop("plot_height must be an unquoted positive integer")
    }} else {
      # ggplots::facet_wrap has argument "space";
      # space = "free_x" is ignored by plotly
      # (known but unsolved issue (08/2022)));
      # the only workaround is to use plotly::subplot;
      # plotly::subplot has no absolute height argument;
      # it can only be given to plotlyOutput;
      # shinybrowser::get_height() can detect height only on server side;
      # ui starts before server and no variable can be passed from server to ui;
      # workaround is the follwing:
      # get standard dev.size and multiply by constant factor;
      # dev.size in RStudio can be manipulated manually and therefore is not
      # reliable
      plot_height <-
        system("R -e 'dev.size(units = c(\"px\"))[2]'", intern = TRUE) %>%
        .[grep("\\[1\\]", .)] %>%
        {as.numeric(base::strsplit(., " ")[[1]][2]) * 1.6}
    }

  ui <- fillPage(
    tags$head(
      tags$style(HTML("
       .shiny-output-error-validation {
         color: #ff0000;
         font-weight: bold;
       }
       ")),
      # https://stackoverflow.com/questions/17325521/r-shiny-display-loading-message-while-function-is-running
      tags$style(type = "text/css", "
       #loadmessage {
         position: fixed;
         top: 0px;
         left: 0px;
         width: 100%;
         padding: 5px 0px 5px 0px;
         text-align: center;
         font-weight: bold;
         font-size: 100%;
         color: #000000;
         background-color: #CCFF66;
         z-index: 105;
       }
       ")),
    shinybrowser::detect(),
    shinyjs::useShinyjs(),
    sidebarPanel(
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                       tags$div("Loading...", id = "loadmessage")),
      awesomeRadio(
        inputId = "rank",
        label = "Rank:",
        choices = colnames(obj@tax_table),
        selected = colnames(obj@tax_table)[1],
        status = "warning"),
      awesomeRadio(
        inputId = "color",
        label = "Aggregate:",
        choices = colnames(obj@sam_data),
        selected = colnames(obj@sam_data)[1],
        status = "warning"),
      awesomeRadio(
        inputId = "facet",
        label = "Facet:",
        choices = c("none", colnames(obj@sam_data)),
        selected = "none",
        status = "warning"),
      hr(),
      sliderInput("rotate", "Rotate x-Axis Text:",
                  min = -90, max = 90, step = 10, value = 0),
      width = 2),
    mainPanel(
      width = 10,
      plotlyOutput("taxa", height = plot_height)
  ))

  server <- function(input, output) {
    output$taxa <- renderPlotly({
      validate(need(input$color != input$facet ,
      "Error: aggregate and facet must not be the same value.
                 Please revise your choice."))
      if(input$facet == "none") facet <- NULL
      else facet <- input$facet
      p <- taxa.plot(obj, input$rank, input$color, facet)
      p[["x"]][["layout"]][["yaxis"]][["range"]] <- c(0,1)
      n <- names(p[["x"]][["layout"]]) %>% stringr::str_detect('xaxis') %>%
        names(p[["x"]][["layout"]])[.]
      for(i in n) {
        p[["x"]][["layout"]][[i]][["tickangle"]] <- as.numeric(input$rotate)
      }
      p <- FacileViz::unify_legend(p)
      p
    })
  }
  shinyApp(ui, server)
}


taxa.abund <- function(obj, tax, add_meta = TRUE){
  if(!obj@otu_table@taxa_are_rows) x <- t(as(otu_table(obj), "matrix"))
  else x <- as(otu_table(obj), "matrix")
  x <- aggregate(x, data.frame(tax_table(obj)[ , tax, drop = FALSE]), sum)
  if(add_meta == TRUE){
    x <- t(x) %>% as.data.frame() %>% `colnames<-`(.[1,]) %>%
      {. <- .[-1,]} %>% lapply(., as.numeric)
    x <- cbind(x, sample_data(obj))
  }
  return(x)
}


taxa.plot <- function (obj, tax, variable, facet_var = NULL, add_meta = TRUE){
  data <- taxa.abund(obj = obj, tax = tax, add_meta = add_meta)
  d <- data[, c(1:(ncol(data) - ncol(sample_data(obj))),
                which(colnames(data) == variable),
                which(colnames(data) == facet_var))]
  d <- reshape2::melt(d, id.vars = c(variable, facet_var))
  f <- function(...) as.list(substitute(...()))
  arglist <- f(data = d, formula = paste(variable, "~ ."), fun.aggregate = sum)
  arglist$formula <- eval(arglist$formula)
  x <- do.call("recast", args = arglist)
  colnames(x) <- c(variable, "total")
  z <- merge(d, x)
  pal <- c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  times <- ceiling(length(unique(z$variable))/9)

  if (is.null(facet_var)) {
    colnames(z) <- c("data", "variable", "value", "total")
    z <- z %>% group_by(data, variable) %>% mutate(rel_ab = sum(value/total*100))
    p <- z %>%
      ggplot(., aes(x = data, y = value/total, fill = variable, col = variable,
                    text = paste(data, "\n",
                                 "Rank:                    ", variable, "\n",
                                 "Rel. Abundance:   ", rel_ab, "%"))) +
      geom_col() + theme_bw() + xlab("") + ylab("") +
      scale_fill_manual(values = rep(pal, times)) +
      scale_color_manual(values = rep(pal, times))
    p <- ggplotly(p, tooltip = "text")
  }

  if (!is.null(facet_var)) {
    colnames(z) <- c("data", "facet", "variable", "value", "total")
    z <- z %>% group_by(facet, data) %>% mutate(facet_total = sum(value))
    z <- z %>% group_by(facet) %>% mutate(plot_width = length(unique(data)))
    plot_widths <- z %>% group_by(facet) %>% summarise(plot_width)
    plot_widths <- unique(as.data.frame(plot_widths))[,2]
    plot_widths <- plot_widths/sum(plot_widths)
    z <- z %>% group_by(facet,data,variable) %>%
      mutate(rel_ab = sum(value/facet_total*100))
    p <- z %>%
      split(.$facet) %>%
      map(function(x) {
        (ggplot(x, aes(x = data, y = value/facet_total,
                       fill = variable, col = variable,
                       text = paste(data, "\n",
                                    "Rank:                    ", variable, "\n",
                                    "Rel. Abundance:   ", rel_ab, "%"))) +
           geom_col() +
           scale_fill_manual(values = rep(pal, times)) +
           scale_color_manual(values = rep(pal, times)) +
           theme_bw() +
           facet_grid(~facet, scales = "free_x", space = "free_x")) %>%
          ggplotly(tooltip = "text")
      }) %>%
      subplot(margin = 0.005, shareY = T, widths = plot_widths)
  }
  return(p)
}
