facet_layout<-function (plot = NULL, facets = NULL, nrow = 2, ncol = 2, scales = "fixed") 
{
  if (is.null(plot)) {
    stop("Argument \"plot\" required")
  }
  if (is.null(facets)) {
    message("Argument \"facets\" not provided. Ploting single panel")
    return(plot)
  }
  if (!all(facets %in% colnames(plot$data))) {
    stop(paste("The facets:", facets, "could not be found in the data"))
  }
  if (is.null(ncol) | is.null(nrow)) {
    stop("Arguments \"ncol\" and \"nrow\" required")
  }
  n_panel_tot <- nrow(unique(plot$data[, facets, drop = FALSE]))
  n_layout <- ncol * nrow
  if (n_panel_tot > n_layout) {
    stop("nrow * ncol >= n is not TRUE, use \"facet_multiple()\" instead")
  }
  n_missing <- n_layout - n_panel_tot
  nrow_last <- max(which(n_panel_tot >= seq(1, n_layout, by = ncol)))
  panel_last <- n_panel_tot - ncol * (nrow_last - 1)
  if (n_missing == 0 || nrow_last == nrow & nrow != 1) {
    plot <- plot + facet_wrap(facets = facets, ncol = ncol, 
                              scales = scales)
    return(plot)
  }
  plot$data[, "panel_to_drop"] <- interaction(plot$data[, 
                                                        facets])
  levels(plot$data[, "panel_to_drop"]) <- c(levels(plot$data[, 
                                                             "panel_to_drop"]), paste0("empty", 1:n_missing))
  plot <- plot + facet_wrap(facets = "panel_to_drop", ncol = ncol, 
                            drop = TRUE, scales = scales)
  drop_grobs <- function(sep = "", ...) {
    if (scales %in% c("free", "free_x")) {
      tmp <- paste("axis_b", (n_panel_tot + 1):n_layout, 
                   sep = sep)
    }
    else if (panel_last < ncol) {
      tmp <- paste("axis_b", which((1:(nrow * ncol)) > 
                                     (panel_last + (nrow - 1) * ncol)), sep = sep)
    }
    else {
      tmp <- NULL
    }
    tmp2 <- as.vector(outer(c("panel", "strip_t", "axis_l"), 
                            (n_panel_tot + 1):n_layout, paste, sep = sep))
    return(c(tmp, tmp2))
  }
  g <- ggplotGrob(plot)
  g$grobs[names(g$grobs) %in% drop_grobs()] <- NULL
  g$layout <- g$layout[!g$layout$name %in% drop_grobs(sep = "-"), 
                       ]
  g$layout[g$layout$name %in% paste0("axis_b-", (1:panel_last) + 
                                       ((nrow - 1) * ncol)), "b"] <- seq(1, 40, by = 4)[nrow_last + 
                                                                                          1]
  grid::grid.newpage()
  grid::grid.draw(g)
}
