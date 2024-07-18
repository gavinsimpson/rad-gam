# Analyse the serpentine grassland data of Green et al (2009)
#
# Green, J., Harte, J., & Ostling, A. (2019). Data from: Species richness,
# endemism, and abundance patterns: tests of two fractal models in a serpentine
# grassland [Data set]. Zenodo. https://doi.org/10.6078/D1MQ2V
#
# Data set downloaded 26 Feb 2024
#
# Observations are on a 16x16 grid, where the first [row,column], $a_{11}$ of
# the grid is in row 1 of the data, while $a_{12}$ is in the second row of the
# data, etc, hence the data are stored in row-wise, filling the rows
# sequentially.

# packages
pkgs <- c("gratia", "mgcv", "scam", "ggplot2", "MRFtools", "dplyr", "readr",
    "here", "stars", "spdep", "MRFtools", "tidyr", "mboost")
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

spp <- readr::read_csv(here("data/serpentine-grasslands/Serpentine.csv"),
    col_type = "dddddddddddddddddddddddd")
spp_list <- readr::read_tsv(here("data/serpentine-grasslands",
    "Serpentine_SpeciesList.txt"), col_types = "cc",
    col_names = c("code", "name"))

# grid
g <- data.frame(expand.grid(col = seq(0, 3.75, by = 0.25) + (0.25 / 2),
  row = seq(0, 3.75, by = 0.25) + (0.25 / 2)),
z = 1:256) |>
  st_as_stars(xy = c("col", "row"), dims = c("col", "row"))

ggplot() + geom_stars(data = g) +
  scale_y_reverse()

# gets grid in the right order as data
g |> as.data.frame() |> arrange(z)

eps <- sqrt(.Machine$double.eps)

g |> st_as_sf() |>
  poly2nb(queen = TRUE, snap = eps)

pen <- g |>
  st_as_sf() |>
  mrf_penalty(node_labels = paste(rep(1:16, each = 16), rep(1:16, times = 16),
    sep = ","))

as.matrix(pen)[1:17, 1:17]

# create RAD, need to work row-wise and apply rank over the vector

# We pivot spp longer initially, then group_by() row,col combination and then
# rank()
serp_rad <- spp |>
  bind_cols(data.frame(row = rep(1:16, each = 16),
    col = rep(1:16, times = 16))) |>
  relocate(row, col, .before = 1L) |>
  pivot_longer(!c(row, col), names_to = "species", values_to = "abundance") |>
  mutate(id = paste(row, col, sep = ","), f_id = factor(id)) |>
  group_by(id) |>
  mutate(rank = rank(-abundance, ties.method = "first")) |>
  ungroup()

# Can now plot the rad
serp_rad |>
  ggplot(aes(x = rank, y = abundance)) +
  geom_point(size = 0.5) +
  facet_grid(row ~ col)

# GAMs
ctrl <- gam.control(nthreads = 3, trace = TRUE)

# model as spatially smooth but not constraining f(rank) marginal to be
# monotonic decreasing
m_rad1 <- bam(abundance ~ te(col, row, rank, d = c(2, 1), bs = c("ds", "cr"),
  k = c(250, 5)),
data = serp_rad, family = nb(), method = "fREML", control = ctrl,
discrete = TRUE)

# model as spatially smooth via MRF but not constraining f(rank) marginal to be
# monotonic decreasing
m_rad2 <- bam(abundance ~ te(f_id, rank, bs = c("mrf", "cr"),
  k = c(256, 5), xt = list(penalty = pen)),
data = serp_rad, family = nb(), method = "fREML", control = ctrl,
discrete = FALSE)

# scams seems the tensor product constructions are limited
# instead here we fit a factor-by convex decreasing smooth with ranef for plot
# mean abundance.
m_rad3 <- scam(abundance ~ s(f_id, bs = "re") +
  s(rank, by = f_id, bs = "mdcxBy", k = 5),
data = as.data.frame(serp_rad), family = poisson())

# this variant allows for the plot means to vary spatially via an MRF
m_rad4 <- scam(abundance ~ s(f_id, bs = "mrf", xt = list(penalty = pen)) +
  s(rank, by = f_id, bs = "mdcxBy", k = 5),
data = as.data.frame(serp_rad), family = poisson())

# mboost also allows monotonically constrained smooths also, and we might be
# able to combine them into tensor product constructions via %X%
# start with a simple additive model
m_rad5 <- gamboost(abundance ~
  bmono(rank, constraint = "decreasing", knots = 5, df = 6) +
  bbs(col, row, knots = 7, df = 6),
data = serp_rad, family = NBinomial())

# next we fit additive but MRF instead of a f(col,row)
m_rad6 <- gamboost(abundance ~
  bmono(rank, constraint = "decreasing", knots = 5, df = 6) +
  bmrf(f_id, bnd = as.matrix(pen)),
data = serp_rad, family = NBinomial())

# next we have a boosted GAM version of m_rad1
m_rad7 <- gamboost(abundance ~
  bmono(rank, constraint = "decreasing", knots = 5, df = 6) %X%
  bbs(col, row, knots = 7, df = 6),
data = serp_rad, family = NBinomial())

# finally we have a boosted GAM version of m_rad1
m_rad8 <- gamboost(abundance ~
  bmono(rank, constraint = "decreasing", knots = 5, df = 6) %X%
  bmrf(f_id, bnd = as.matrix(pen)),
data = serp_rad, family = NBinomial())
