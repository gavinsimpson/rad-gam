# Loch of the Lowes

# packages
pkgs <- c("gratia", "mgcv", "ggplot2", "MRFtools", "dplyr", "readr",
    "here", "tidyr", "cepreader", "scam")
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

# read the data - it's in CEP format so we need a special reader
lowe2 <- readCEP(here("data/loch-of-lowes/lowe2.cep")) |>
  tibble::rownames_to_column(var = "sample") |>
  as_tibble() |>
  mutate(depth = gsub("^LOWE", "", sample),
    depth = as.numeric(depth)) |>
  relocate(depth, .after = 1L)

lowe2_rad <- lowe2 |>
  pivot_longer(!c(depth, sample),
    names_to = "species", values_to = "abundance") |>
  group_by(sample) |>
  mutate(rank = rank(-abundance, ties.method = "first"),
    f_sample = factor(sample),
    f_depth = factor(depth, levels = lowe2$depth)) |>
  ungroup()

lowe2_rad |>
  ggplot(aes(x = rank, y = abundance, group = f_sample, colour = depth)) +
  geom_point() +
  facet_wrap(~ f_depth) +
  scale_colour_viridis_c(direction = -1)

ctrl <- gam.control(trace = TRUE)

lowe2_m1 <- bam(abundance ~
  s(rank, bs = "cr", k = 6) +
  s(rank, f_depth, bs = "fs", k = 6),
data = lowe2_rad,
family = nb(),
method = "fREML",
discrete = TRUE,
control = ctrl,
nthreads = 3
)

lowe2_m2 <- bam(abundance ~
    s(rank, bs = "cr", k = 6) +
    s(depth, bs = "cr", k = 10) +
    ti(rank, depth, bs = c("cr", "cr"), k = c(10, 10)),
  data = lowe2_rad,
  family = nb(),
  method = "fREML",
  discrete = TRUE,
  control = ctrl,
  nthreads = 3
)

lowe2_m3 <- bam(abundance ~
    s(rank, bs = "ad", k = 10) +
    s(depth, bs = "cr", k = 10) +
    ti(rank, depth, bs = c("cr", "cr"), k = c(10, 10)),
  data = lowe2_rad,
  family = nb(),
  method = "fREML",
  discrete = TRUE,
  control = ctrl,
  nthreads = 3
)

AIC(lowe2_m1, lowe2_m2, lowe2_m3)

lowe_ds_m1 <- data_slice(
  lowe2_m1,
  rank = evenly(rank),
  f_depth = evenly(f_depth)
)

lowe_ds_m2 <- data_slice(
  lowe2_m2,
  rank = evenly(rank),
  depth = unique(depth)
)

fv_lowe_m1 <- fitted_values(lowe2_m1, data = lowe_ds_m1)
fv_lowe_m2 <- fitted_values(lowe2_m2, data = lowe_ds_m2)
fv_lowe_m3 <- fitted_values(lowe2_m3, data = lowe_ds_m2)

lowe2_rad |>
  ggplot(aes(x = rank, y = abundance, group = f_depth)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ f_depth) +
  geom_line(data = fv_lowe_m1,
    aes(y = .fitted, x = rank, group = f_depth),
    colour = "red") +
  labs(x = "Species rank", y = "Abundance")

lowe2_rad |>
  ggplot(aes(x = rank, y = abundance, group = depth)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ depth) +
  geom_line(data = fv_lowe_m2,
    aes(y = .fitted, x = rank, group = depth),
    colour = "red") +
  labs(x = "Species rank", y = "Abundance")

lowe2_rad |>
  ggplot(aes(x = rank, y = abundance, group = depth)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~ depth) +
  geom_line(data = fv_lowe_m3,
    aes(y = .fitted, x = rank, group = depth),
    colour = "red") +
  labs(x = "Species rank", y = "Abundance")

lowe2_rad |>
  filter(depth == 0) |>
  ggplot(aes(x = rank, y = abundance)) +
  geom_point() +
  geom_ribbon(
    data = fv_lowe_m3 |>
      filter(depth == 0),
    aes(
      ymin = .lower_ci, ymax = .upper_ci, y = .fitted
    ),
    fill = "red", alpha = 0.2
  ) +
  geom_line(data = fv_lowe_m3 |> filter(depth == 0),
    aes(y = .fitted, x = rank, group = depth),
    colour = "red") +
  labs(x = "Species rank", y = "Abundance")

draw(lowe2_m3, rug = FALSE)

lowe_ds_rank_ts <- data_slice(
  lowe2_m3,
  rank = evenly(rank, by = 1),
  depth = evenly(depth)
)

fv_lowe_rank_ts <- fitted_values(lowe2_m3, data = lowe_ds_rank_ts)

brks <- seq(2, 10, by = 2)
fv_lowe_rank_ts |>
  filter(rank <= 10) |> # the top 10 species
  ggplot(aes(x = depth, y = .fitted, group = rank)) +
  geom_ribbon(
    aes(
      ymin = .lower_ci,
      ymax = .upper_ci,
      y = .fitted,
      x = depth,
      fill = rank),
    alpha = 0.2) +
  geom_line(
    aes(
      colour = rank
    ),
    linewidth = 1
  ) +
  scale_x_reverse() +
  scale_fill_viridis_c(
    breaks = brks,
    option = "inferno",
    direction = -1
  ) +
  scale_colour_viridis_c(
    breaks = brks,
    option = "inferno",
    direction = -1
  ) +
  labs(y = "Abundance", x = "Sample depth", fill = "Rank", colour = "Rank")

# trying Dave's suggestion to trick scam into doing something
# this doesn't work - it fits but something is broken so
# summary(), predict() etc just fail.
# need to email Natalya about it
lowe2_s1 <- scam(abundance ~
  s(rank, f_depth, bs = "fs", k = 6, xt = list(bs = "mpd")),
data = lowe2_rad,
family = poisson(),
control = scam.control(trace = TRUE)
)

lowe_ds_s1 <- data_slice(
  lowe2_s1,
  rank = evenly(rank),
  f_depth = evenly(f_depth)
)

fv_lowe_s1 <- fitted_values(lowe2_s1, data = lowe_ds_s1)
