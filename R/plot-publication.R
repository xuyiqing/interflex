# Presentation-only styling. Layer data, mappings, model estimates and confidence
# limits remain untouched. Clone mutable ggproto layers before changing widths.
.interflex_plot_preset <- function(preset) {
    if (is.null(preset)) return(FALSE)
    if (!is.character(preset) || length(preset) != 1L ||
        is.na(preset) || preset != "publication") {
        stop('"preset" must be NULL or "publication".')
    }
    TRUE
}

.interflex_publication_plot <- function(p, axis.size = 18, tick.size = 12,
                                       add.zero = TRUE) {
    if (!inherits(p, "ggplot")) return(p)
    layers <- lapply(p$layers, function(layer) {
        copy <- ggplot2::ggproto(NULL, layer)
        if (inherits(copy$geom, "GeomErrorbar")) {
            copy$aes_params$width <- 0
            copy$geom_params$width <- 0
        }
        copy
    })
    if (add.zero) {
        # Replace only the default horizontal zero reference. Other user-added
        # reference layers and all estimate/interval data are preserved.
        is.zero <- vapply(layers, function(layer) {
            inherits(layer$geom, "GeomHline") &&
                is.data.frame(layer$data) &&
                "yintercept" %in% names(layer$data) &&
                length(layer$data$yintercept) > 0L &&
                all(is.finite(layer$data$yintercept) & layer$data$yintercept == 0)
        }, logical(1))
        layers <- c(list(ggplot2::geom_hline(yintercept = 0, colour = "black",
                                            linewidth = 0.65)), layers[!is.zero])
    }
    p$layers <- layers
    p <- p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
        ggplot2::theme_bw(base_size = tick.size) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = axis.size, colour = "black"),
            axis.text = ggplot2::element_text(size = tick.size, colour = "black"),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 9)),
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "gray45", linewidth = 0.35),
            plot.margin = ggplot2::margin(14, 18, 14, 14),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = tick.size),
            legend.key.size = grid::unit(0.9, "lines"),
            legend.background = ggplot2::element_rect(fill = "white", colour = NA),
            legend.margin = ggplot2::margin(2, 3, 2, 3)
        )
    if (utils::packageVersion("ggplot2") >= "3.5.0") {
        p <- p + ggplot2::theme(legend.position = "inside",
                                legend.position.inside = c(0.02, 0.98),
                                legend.justification = c(0, 1))
    } else {
        p <- p + ggplot2::theme(legend.position = c(0.02, 0.98),
                                legend.justification = c(0, 1))
    }
    p
}
