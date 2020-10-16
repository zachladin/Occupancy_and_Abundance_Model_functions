facet_multiple<- function(plot = NULL, facets = NULL, ncol = 2, nrow = 2, page = NULL, 
          scales = "fixed") 
{
  require(ggplot2)
  
  if (is.null(plot)) {
    stop("Argument \"plot\" required")
  }
  if (is.null(facets)) {
    message("Argument \"facets\" not provided. Ploting single panel")
    return(plot)
  }
  if (!all(facets %in% colnames(plot$data))) {
    stop("The facet variable ", paste(facets, collapse = ", "), 
         " not found")
  }
  if (is.null(ncol) | is.null(nrow)) {
    stop("Arguments \"ncol\" and \"nrow\" required")
  }
  n_panel_tot <- nrow(unique(plot$data[, facets, drop = FALSE]))
  n_layout <- ncol * nrow
  n_pages <- ceiling(n_panel_tot/n_layout)
  if (!is.null(page) && page > n_pages) {
    stop("Argument \"page\" > last page (p=", n_pages, ")")
  }
  plot <- plot + facet_wrap(facets = facets, ncol = ncol, 
                            scales = scales)
  if (n_pages == 1) {
    return(plot)
  }
  data <- plot$data
  title <- plot$labels$title
  if (!scales %in% c("free", "free_x") && is.numeric(eval(plot$mapping$x, 
                                                          data)) && length(grep("xmax", plot$scales$scales, fixed = TRUE)) == 
      0) {
    plot$coordinates$limits$x <- range(eval(plot$mapping$x, 
                                            data))
  }
  if (!scales %in% c("free", "free_y") && is.numeric(eval(plot$mapping$y, 
                                                          data)) && length(grep("ymax", plot$scales$scales, fixed = TRUE)) == 
      0) {
    plot$coordinates$limits$y <- range(eval(plot$mapping$y, 
                                            data))
  }
  data$groups <- findInterval(unclass(interaction(data[, facets])), 
                              seq(from = 1, by = n_layout, length.out = n_pages)[-1]) + 
    1
  if (is.null(page)) {
    draw_page <- seq_along(1:n_pages)
  }
  else {
    draw_page <- page
  }
  for (i in draw_page) {
    plot <- plot %+% data[data$groups == i, ] + ggtitle(label = bquote(atop(bold(.(title)), 
                                                                            atop(italic(Page ~ .(i) ~ of ~ .(n_pages))))))
    if (i == n_pages) {
      plot <- facet_layout(plot = plot, facets = facets, 
                           ncol = ncol, nrow = nrow, scales = scales)
    }
    if (!is.null(plot)) {
      print(plot)
    }
  }
}