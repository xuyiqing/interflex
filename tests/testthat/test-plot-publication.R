publication_fit <- local({
    cached <- NULL
    function() {
        if (!is.null(cached)) return(cached)
        set.seed(341)
        dat <- expand.grid(X = 0:2, D = 0:1, replication = seq_len(40))
        dat$Y <- dat$D * (0.5 + dat$X / 3) + rnorm(nrow(dat))
        cached <<- suppressWarnings(suppressMessages(interflex(
            data = dat, Y = "Y", D = "D", X = "X", estimator = "linear",
            gate = TRUE, vartype = "delta", nboots = 5, figure = FALSE
        )))
        cached
    }
})

test_that("publication style preserves estimates and CI endpoints", {
    fit <- publication_fit()
    before <- serialize(fit$g.est, NULL)
    ordinary <- plot(fit, show.all = TRUE, Xdistr = "none")[[1]]
    p <- plot(fit, preset = "publication", show.all = TRUE,
              Xdistr = "none", main = "Put this in the manuscript")[[1]]
    expect_identical(serialize(fit$g.est, NULL), before)
    expect_null(p$labels$title)
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
    expect_identical(p$labels$y, "Estimate")
    expect_equal(p$theme$axis.title$size, 18)
    expect_equal(p$theme$axis.text$size, 12)

    bars <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomErrorbar"), TRUE))
    old.bars <- which(vapply(ordinary$layers, function(l) inherits(l$geom, "GeomErrorbar"), TRUE))
    expect_gt(length(bars), 0)
    built <- ggplot2::ggplot_build(p)$data
    old.built <- ggplot2::ggplot_build(ordinary)$data
    for (i in seq_along(bars)) {
        b <- built[[bars[i]]]
        expect_equal(b$xmin, b$xmax)
        expect_equal(b$ymin, old.built[[old.bars[i]]]$ymin)
        expect_equal(b$ymax, old.built[[old.bars[i]]]$ymax)
    }
    # Styling a new plot must not mutate a previously returned plot's layers.
    expect_true(any(old.built[[old.bars[1]]]$xmin != old.built[[old.bars[1]]]$xmax))
    zeros <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomHline"), TRUE))
    expect_length(zeros, 1)
    expect_equal(p$layers[[zeros]]$aes_params$colour, "black")
})

test_that("publication keeps caller labels, colors, sizes and limits", {
    p <- plot(publication_fit(), preset = "publication", show.all = TRUE,
              Xdistr = "none", xlab = "Custom X", ylab = "Custom Y",
              line.color = "purple", cex.lab = 1.4, cex.axis = 1,
              xlim = c(0, 2), ylim = c(-1, 2))[[1]]
    expect_identical(p$labels$x, "Custom X")
    expect_identical(p$labels$y, "Custom Y")
    expect_equal(p$theme$axis.title.x$size, 21)
    expect_equal(p$theme$axis.text$size, 15)
    expect_equal(p$coordinates$limits$x, c(-0.08, 2.08))
    expect_equal(p$coordinates$limits$y, c(-1, 2))
    points <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE))
    expect_true(all(vapply(p$layers[points], function(l) identical(l$aes_params$fill, "purple"), TRUE)))
})

test_that("pooled publication renderer and file output work", {
    fit <- publication_fit()
    file <- tempfile(fileext = ".png")
    on.exit(unlink(file))
    p <- suppressMessages(plot(fit, pool = TRUE, preset = "publication",
                               Xdistr = "none", file = file, width = 6, height = 4))
    expect_s3_class(p, "ggplot")
    expect_true(file.exists(file))
    expect_gt(file.info(file)$size, 0)
    g <- ggplot2::ggplotGrob(p)
    inside.name <- if (utils::packageVersion("ggplot2") >= "3.5.0") "guide-box-inside" else "guide-box"
    inside <- which(g$layout$name == inside.name)
    expect_length(inside, 1)
    expect_false(inherits(g$grobs[[inside]], "zeroGrob"))
    # Without the legacy pooled legend carrier, inspect the styled plot directly.
    raw <- suppressMessages(plot(fit, pool = TRUE, preset = "publication",
                                 Xdistr = "none", show.subtitles = FALSE))
    expect_null(raw$labels$title)
    expect_equal(raw$theme$axis.title$size, 18)
    expect_identical(raw$labels$y, "Estimate")
    # Default palette is grayscale; an explicit palette remains available.
    lines <- which(vapply(raw$layers, function(l) inherits(l$geom, "GeomLine"), TRUE))
    built <- ggplot2::ggplot_build(raw)$data
    expect_true(all(built[[lines[1]]]$colour == "#000000"))
    custom <- suppressMessages(plot(fit, pool = TRUE, preset = "publication",
                                    Xdistr = "none", show.subtitles = FALSE,
                                    color = c("gray70", "purple")))
    custom.lines <- which(vapply(custom$layers, function(l) inherits(l$geom, "GeomLine"), TRUE))
    expect_true(all(ggplot2::ggplot_build(custom)$data[[custom.lines[1]]]$colour == "purple"))
})

test_that("NULL preset preserves the default and invalid presets fail early", {
    fit <- publication_fit()
    p <- plot(fit, show.all = TRUE, Xdistr = "none")
    q <- plot(fit, show.all = TRUE, Xdistr = "none", preset = NULL)
    expect_equal(ggplot2::ggplot_build(p[[1]])$data,
                 ggplot2::ggplot_build(q[[1]])$data)
    expect_error(plot(fit, preset = "unknown"), "preset")
    expect_error(plot(fit, preset = NA_character_), "preset")
})
