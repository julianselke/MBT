#' Plot QC Summaries
#'
#' Various plotting options for results of readQC.
#'
#' @md
#' @param data A list returned by readQC().
#' @param option One of:
#' * "qp": Quantile plot
#' * "mq": Mean quality per cycle.
#' * "rc": Read count per sample.
#' * "bp": Quantile boxplot (similar to quality plot in qzv artifact produced
#' by **qiime demux summarize**)
#' * "hm": Quality heatmap
#' @param meta_data_col Character specifying the the metadata column to color
#' read counts by.
#' @param fill.log A boolean. Scale fill of option "hm".
#' @param add.mean A boolean. Plot means on top of heatmap ("hm")
#' @param q Integer. One of 1:length(qs)
#' @param qs Numeric vector. Specifies quantiles calculated in readQC().
#' @param ... Passed to \link[viridis]{scale_fill_viridis_c}
#' @return A ggplot object.
#' @export
#' @keywords MBT
#' @examples plotQC(data)
#' @import ggplot2
#' @import dplyr
#' @import forcats
#' @import ggsci



plotQC <- function(data,
                   option,
                   meta_data_col,
                   fill.log = FALSE,
                   add.mean = TRUE,
                   q = 3,
                   qs = c(0.09, 0.25, 0.5, 0.75, 0.91), ...) {
  switch(option,
         "qp" = {
           .dat <- data[["read.quants"]] %>%
             .[.$variable == as.character(qs[q]), ] %>%
             .[order(.[, 1], .[, 2]), ]
           ggplot(data = .dat) +
             theme_bw() +
             geom_line(aes(x = Cycle, y = value, group = SampleID),
                       alpha = 0.1) +
             scale_x_continuous(breaks = seq(20, max(.dat$Cycle), 20)) +
             scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
             coord_cartesian(expand = FALSE) +
             theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title = element_blank(),
               legend.position = "none"
             )
         },
         "mq" = {
           ggplot(data[["read.means"]]) +
             geom_line(aes(group = Var1, x = Var2, y = value), alpha = 0.1) +
             theme_bw() +
             scale_x_continuous(breaks = seq(20, max(data[[2]]$Var2), 20)) +
             scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
             coord_cartesian(expand = FALSE) +
             theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title = element_blank(),
               legend.position = "none"
             )
         },
         "rc" = {
           if (!meta_data_col %in% colnames(data[["read.counts"]])) {
             did_you_mean(colnames(data[["read.counts"]]), meta_data_col)
           }
           data[["read.counts"]] %>%
             mutate(name = fct_reorder(Sample, Read.Counts, .desc = TRUE)) %>%
             ggplot() +
             geom_col(aes(
               x = name, y = Read.Counts,
               fill = get(meta_data_col)
             )) +
             theme_bw() +
             scale_y_continuous(
               expand = expansion(mult = c(0, 0.1)),
               breaks = seq(0, 1000000, 10000)
             ) +
             theme(
               axis.text.x = element_text(angle = -90),
               legend.position = "bottom"
             ) +
             ylab("Read Count") +
             xlab("Sample") +
             {
               if (!is.numeric(data[["read.counts"]][[eval(meta_data_col)]])) {
                 scale_fill_d3(name = eval(meta_data_col))
               } else {
                 scale_fill_viridis_c(name = eval(meta_data_col), ...)
               }
             }
         },
         "bp" = {
           .dat <- data[["read.quants"]]
           ggplot() +
             geom_boxplot(
               aes(
                 x = .dat[.dat$variable == "0.09", ]$Cycle %>%
                   unique() %>%
                   sort(),
                 group = .dat[.dat$variable == "0.09", ]$Cycle %>%
                   unique() %>%
                   sort(),
                 ymin = .dat[.dat$variable == "0.09", ] %>%
                   group_by(Cycle) %>%
                   summarize(Mean = mean(value, na.rm = TRUE)) %>%
                   .$Mean,
                 lower = .dat[.dat$variable == "0.25", ] %>%
                   group_by(Cycle) %>%
                   summarize(Mean = mean(value, na.rm = TRUE)) %>%
                   .$Mean,
                 middle = .dat[.dat$variable == "0.5", ] %>%
                   group_by(Cycle) %>%
                   summarize(Mean = mean(value, na.rm = TRUE)) %>%
                   .$Mean,
                 upper = .dat[.dat$variable == "0.75", ] %>%
                   group_by(Cycle) %>%
                   summarize(Mean = mean(value, na.rm = TRUE)) %>%
                   .$Mean,
                 ymax = .dat[.dat$variable == "0.91", ] %>%
                   group_by(Cycle) %>%
                   summarize(Mean = mean(value, na.rm = TRUE)) %>%
                   .$Mean
               ),
               stat = "identity",
               fill = "steelblue",
               linetype = 1,
               color = "gray20",
               lwd = 0.2
             ) +
             theme_bw() +
             scale_x_continuous(breaks = seq(20, 1000, 20)) +
             scale_y_continuous(
               breaks = seq(0, 40, 5),
               limits = c(0, 40),
               minor_breaks = seq(1, 40, 1)
             ) +
             coord_cartesian(expand = FALSE) +
             theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "none"
             )
         },
         "hm" = {
           ggplot() +
             {
               if (fill.log) {
                 geom_tile(
                   data = data[["quality.tiles"]],
                   aes(x = Cycle, y = Score, fill = log(Count))
                 )
               }
             } +
             {
               if (!fill.log) {
                 geom_tile(
                   data = data[["quality.tiles"]],
                   aes(x = Cycle, y = Score, fill = Count)
                 )
               }
             } +
             {
               if (add.mean) {
                 geom_line(
                   data = data[["read.means"]],
                   aes(group = Var1, x = Var2, y = value),
                   color = "red", alpha = 0.1
                 )
               }
             } +
             scale_fill_viridis_c(
               direction = -1,
               option = "mako",
               na.value = "white"
             ) +
             theme_bw() +
             scale_x_continuous(breaks = seq(20, 1000, 20)) +
             scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
             coord_cartesian(expand = FALSE) +
             theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.position = "none"
             )
         }
  )
}
