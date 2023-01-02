#' Rarefaction Overview
#'
#' Interactively assess effects of rarefaction depth on features and samples.
#'
#' @md
#' @param obj A _phyloseq_ object.
#' @return Starts a Shiny App.
#' @keywords MBT
#' @export
#' @examples
#' rarefySH(GlobalPatterns)
#' @import ggplot2
#' @import dplyr
#' @import vegan
#' @import shiny
#' @import phyloseq

rarefySH <- function(obj) {
  # taxa rows?
  ft <- otu_table(obj) %>%
    t() %>%
    data.frame() %>%
    `class<-`(., "data.frame")
  md <- sample_data(obj) %>%
    data.frame() %>%
    `class<-`(., "data.frame")
  raremax <- max(rowSums(ft))
  message("Computing Rarefaction Curve...")
  rmq <<- vegan::rarecurve(ft,
                           step = raremax / 200,
                           sample = raremax,
                           tidy = TRUE) %>%
    group_by(Site) %>%
    mutate(Max = max(Sample)) %>%
    data.frame()
  m <<- rmq[, c(1, 4)] %>%
    unique.data.frame() %>%
    `colnames<-`(c(colnames(md)[1], "max")) %>%
    merge(., md, colnames(md)[1])
  pal <- c(
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
  )
  if (interactive()) {
    ui <- fluidPage(
      fluidRow(plotOutput("depthPlot")),
      fluidRow(column(12,
                      align = "center",
                      sliderInput("depth", "Rarefaction Depth:",
                                  min = 0,
                                  max = max(rmq$Sample), value = raremax,
                                  width = "50%"
                      )
      )),
      fluidRow(column(12,
                      align = "center",
                      radioButtons(
                        inputId = "group", inline = TRUE,
                        label = "Select Metadata Variable:",
                        choices = colnames(md[3:ncol(md)])
                      )
      )),
      fluidRow(column(12, align = "center", plotOutput("rarePlot")))
    )

    server <- function(input, output) {
      output$depthPlot <- renderPlot({
        ggplot() +
          geom_line(
            data = rmq,
            aes(x = Sample, y = Species, group = Site),
            color = "gray"
          ) +
          theme_bw() +
          theme(legend.position = "none") +
          geom_line(
            data = rmq[rmq$Max >= input$depth, ],
            aes(x = Sample, y = Species, group = Site),
            color = "#636EFA"
          ) +
          geom_vline(xintercept = input$depth)
      })
      output$rarePlot <- renderPlot({
        p <- ggplot(m) +
          geom_col(aes(x = m[[input$group]], y = max), fill = "gray") +
          geom_col(
            data = m[m$max >= input$depth, ],
            aes(
              x = m[m$max >= input$depth, ][[input$group]],
              y = max,
              fill = m[m$max >= input$depth, ][[input$group]]
            )
          ) +
          geom_hline(yintercept = input$depth) +
          theme_bw() +
          theme(
            legend.position = "none",
            axis.title = element_blank(),
            plot.margin = unit(c(0, 0, 1, 0), "cm")
          )
        if (length(levels(as.factor(m[[input$group]]))) <= 10) {
          p + scale_fill_manual(values = pal)
        } else {
          p
        }
      })
    }
    shinyApp(ui, server)
  }
}
