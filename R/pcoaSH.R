#' PCoA Overview
#'
#' View PCoA plots in an interactive shiny app. Select distance/dissimilarity
#' metrices and subset data.
#'
#' @param obj A _phyloseq_ object.
#' @return Starts a Shiny App.
#' @keywords MBT
#' @export
#' @examples
#' pcoaSH(GlobalPatterns)
#' @import ggplot2 dplyr vegan shiny shinybrowser shinyjs shinyWidgets philr plotly ape phyloseq

pcoaSH <- function(obj){

  distList <- rbind(c("bray","Bray-Curtis"),
                    c("jaccard","Jaccard"),
                    c("u.unifrac","Unweighted UniFrac"),
                    c("w.unifrac","Weighted UniFrac"),
                    c("aitchison","Aitchison"),
                    c("philr","PhILR")) %>%
    as.data.frame() %>%
    `colnames<-`(c("id","name"))

  ui <- fluidPage(
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
    useShinyjs(),
    sidebarPanel(
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                       tags$div("Loading...", id = "loadmessage")),
      awesomeRadio(
        inputId = "dist",
        label = "Metric:",
        choices = distList$name,
        selected = "Bray-Curtis",
        status = "warning"),
      conditionalPanel(
        condition = "input.dist == 'PhILR'",
        pickerInput(
          inputId = "part.weights",
          label = "part.weights",
          choices = c("uniform", "gm.counts", "anorm", "anorm.x.gm.counts",
                      "enorm", "enorm.x.gm.counts")
        ),
        pickerInput(
          inputId = "ilr.weights",
          label = "ilr.weights",
          choices = c("uniform", "blw", "blw.sqrt", "mean.descendants")
        )
      ),
      awesomeRadio(
        inputId = "color",
        label = "Color:",
        choices = sample_variables(obj),
        selected = sample_variables(obj)[1],
        status = "warning"),
      width = 2),
    mainPanel(
      width=8,
      plotlyOutput("pcoa")),
    sidebarPanel(
      h3("Subset Data"),
      width = 2,
      uiOutput("column"),
      uiOutput("level"))
  )

  server <- function(input, output) {

    output$column <- renderUI({
      selectInput(inputId = "column",
                  label = h5("Choose Column"),
                  choices = sample_variables(obj),
                  selected = 1)
    })

    output$level <- renderUI({
      pickerInput(
        inputId = "level",
        label = h5("Choose Level"),
        choices = select(sample_data(obj) %>% `class<-`("data.frame"),
                         input$column) %>% as.vector() %>% .[[1]] %>%
          as.character() %>% unique(),
        selected = select(sample_data(obj) %>% `class<-`("data.frame"),
                          input$column) %>% as.vector() %>% .[[1]] %>%
          as.character() %>% unique(),
        options = list(
          `actions-box` = TRUE),
        multiple = TRUE
      )
    })

    obj_ <- reactive({return(obj)})

    subObj <- reactive({
      .col_ <<- as.symbol(input$column)
      .vals_ <<- input$level
      x <- subset_samples(obj_(), eval(.col_) %in% .vals_)
      return(x)
    })

    idata <- reactive({
      MBT::dist(obj = subObj(), dis = distList$id[distList$name == input$dist],
              pseudocount = 1e-12, ilr.weights = input$ilr.weights,
              part.weights = input$part.weights)
    })

    output$pcoa <- renderPlotly({
      pcoa <- ape::pcoa(idata())
      data <- data.frame(Sample = as.character(rownames(pcoa$vectors)),
                         X = pcoa$vectors[ , 1],
                         Y = pcoa$vectors[ , 2],
                         Z = pcoa$vectors[ , 3])
      data$Sample <- rownames(data)
      colnames(data) <- c("sampleid", "X", "Y", "Z")
      data <- merge(data, sample_data(subObj()), by="row.names")
      axx <- list(title = paste("PCo 1 - ", round(
        pcoa$values$Relative_eig[[1]] * 100, 2), "%", sep = ""),
        tickfont = list(color = c('#ababab')), color = "#777777")
      axy <- list(title = paste("PCo 2 - ", round(
        pcoa$values$Relative_eig[[2]] * 100, 2), "%", sep = ""),
        tickfont = list(color = c('#ababab')), color = "#777777")
      axz <- list(title = paste("PCo 3 - ", round(
        pcoa$values$Relative_eig[[3]] * 100, 2), "%", sep = ""),
        tickfont = list(color = c('#ababab')), color = "#777777")
      plot_ly(type = 'scatter3d', mode = 'markers', data = data, x = ~X,
              y = ~Y, z = ~Z, color = ~get(input$color),
              height = shinybrowser::get_height(),
              colors = c("#00468B","#ED0000","#42B540","#0099B4","#925E9F"),
              marker = list(size = 5, opacity = 1.0), showlegend = T) %>%
        layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz),
               title = "")
    })
  }
  shinyApp(ui, server)
}
