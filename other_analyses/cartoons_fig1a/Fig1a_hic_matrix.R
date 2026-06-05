## Simulate Hi-C data to create a mock Hi-C interaction frequency matrix for Figure 1a. This isn't
## used in any analyses and solely serves as cartoon illustrating Hi-C data

## Define functions --------------------------------------------------------------------------------

# function to simulate a Hi-C interaction matrix with a specific number of bins
simulate_hic <- function(n = 50, seed = 42) {
  
  stopifnot(is.finite(n), n >= 5)   # need a few bins for TADs to make sense
  
  if (!is.null(seed)) set.seed(seed)
  mat <- matrix(0, n, n)
  
  # simulate sistance-dependent decay so that contact frequency falls off as a power law
  for (i in 1:n) {
    for (j in 1:n) {
      d <- abs(i - j)
      mat[i, j] <- 1 / (d + 1)^0.9
    }
  }
  
  # simulate TADs
  n_tads     <- 5
  tad_bounds <- unique(round(seq(1, n, length.out = n_tads + 1)))
  for (k in 1:(length(tad_bounds) - 1)) {
    s <- tad_bounds[k]
    e <- tad_bounds[k + 1]
    mat[s:e, s:e] <- mat[s:e, s:e] * 2.2
  }
  
  # A/B compartments (large-scale plaid pattern)
  comp_width  <- max(1, round(n / 4))
  compartment <- ifelse((1:n %/% comp_width) %% 2 == 0, 1, -1)
  comp_boost  <- outer(compartment, compartment, "*")
  mat <- mat * (1 + 0.25 * comp_boost)
  
  # chromatin loops (bright foci off the diagonal)
  corners <- tad_bounds[-c(1, length(tad_bounds))]
  loops   <- list()
  if (length(corners) >= 2) {
    for (k in 1:(length(corners) - 1)) {
      loops[[length(loops) + 1]] <- c(corners[k], corners[k + 1])
    }
    loops[[length(loops) + 1]] <- c(corners[1], corners[length(corners)])
  }
  for (lp in loops) {
    i <- lp[1]; j <- lp[2]
    if (i != j) {
      mat[i, j] <- mat[i, j] + 0.8
      mat[j, i] <- mat[j, i] + 0.8
    }
  }
  
  # poisson-like noise and enforce symmetry
  counts <- matrix(rpois(n * n, lambda = mat * 50), n, n)
  counts[lower.tri(counts)] <- t(counts)[lower.tri(counts)]
  
  # log-transform for visualization
  logmat <- log1p(counts)
  
  return(logmat)
  
}

# function to plot a Hi-C matrix
plot_hic_matrix <- function(hic, hic_colors, output_file) {
  
  # get number of bins
  n <- nrow(hic)
  
  # print Hi-C matrix to file (.pdf)
  pdf(output_file, width = 8, height = 7.5)
  par(mar = c(4, 4, 3, 6))
  image(1:n, 1:n, hic[, n:1], col = hic_colors,
        xlab = "Genomic bin", ylab = "Genomic bin",
        main = sprintf("Simulated Hi-C Contact Matrix (%d x %d bins)", n, n),
        axes = FALSE, useRaster = FALSE)
  axis(1)
  tick_step <- max(1, round(n / 5)); ticks <- seq(0, n, by = tick_step)
  axis(2, at = ticks, labels = rev(ticks)); box()
  dev.off()

}

## Create mock Hi-C matrix -------------------------------------------------------------------------

# Hi-C matrix with 25 genomic bins
hic25 <- simulate_hic(n = 25)

# color gradient for plotting matrix
hic_colors25 <- colorRampPalette(c(rep("white", 6), "#FCBBA1", "#FC4E2A", "#BD0026", "#8c011d"),
                                 bias = 2.2)(256)

# make plot and save to pdf file
plot_hic_matrix(hic25, hic_colors25, "other_analyses/cartoons_fig1a/fig1a_hic_25bins.pdf")
