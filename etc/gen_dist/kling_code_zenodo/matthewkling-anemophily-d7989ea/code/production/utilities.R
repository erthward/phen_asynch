


# modified version of ggsave
ggs <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
                 scale = 1, width = NA, height = NA, 
                 units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, 
                 add, # new argument
                 ...) 
{
        require(ggplot2)
        require(grid)
        source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/save.r")
        source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/utilities.r")
        dpi <- parse_dpi(dpi)
        dev <- plot_dev(device, filename, dpi = dpi)
        dim <- plot_dim(c(width, height), scale = scale, units = units, 
                        limitsize = limitsize)
        if (!is.null(path)) {
                filename <- file.path(path, filename)
        }
        old_dev <- grDevices::dev.cur()
        dev(filename = filename, width = dim[1], height = dim[2], 
            ...)
        on.exit(utils::capture.output({
                grDevices::dev.off()
                if (old_dev > 1) grDevices::dev.set(old_dev)
        }))
        grid.draw(plot)
        lapply(add, function(x) x)
        invisible()
}
