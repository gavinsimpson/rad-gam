# analysis of the portal rodent data

# packages
pkgs <- c("gratia", "mgcv", "ggplot2", "MRFtools", "dplyr", "readr",
    "here", "MRFtools", "tidyr", "portalr")
vapply(pkgs, library, logical(1L), logical.return = TRUE, character.only = TRUE)

# download the portal data
portal_path <- here("data/portal/")
# download_observations(path = portal_path, timeout = 360) # only needed once

# load data
# tbls <- load_rodent_data(portal_path) # don't need this for RACs
# load rodent abundance data
rodents <- abundance(portal_path, time = "date", level = "plot") |>
  as_tibble()

## is census date the right time variable?
rodents |>
  filter(treatment == "control") |>
  group_by(censusdate) |>
  count() |>
  ungroup() |>
  ggplot(aes(x = censusdate, y = n)) +
  geom_point()
## seems to be OK; but what is the difference in days between observations...
rodents |>
  filter(plot %in% c(4, 11, 14, 17)) |>
  mutate(year = format(censusdate, "%Y")) |>
  group_by(plot, year) |>
  reframe(diff = as.numeric(diff(censusdate))) |>
  ggplot(aes(x = diff, colour = as.factor(plot), group = plot)) +
  geom_density() +
  facet_wrap(~ year)

## OK, seems alright just now, compute the RAC
rodent_rac <- rodents |>
  filter(plot %in% c(4, 11, 14, 17)) |>
  pivot_longer(!c(censusdate, treatment, plot),
    names_to = "species", values_to = "abundance") |>
  group_by(plot, censusdate) |>
  mutate(rank = rank(-abundance, ties.method = "first"),
  year = format(censusdate, "%Y") |> as.numeric(),
  month = format(censusdate, "%m") |> as.numeric(),
  doy = format(censusdate, "%j") |> as.numeric(),
  census_f = factor(censusdate),
  plot_f = factor(plot),
  year_f = factor(year)) |>
  ungroup()

# Can now plot the rac
rodent_rac |>
  ggplot(aes(x = rank, y = abundance,
    colour = censusdate, group = censusdate)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~ plot) +
  scale_colour_viridis_c(option = "magma")

# model
ctrl <- gam.control(trace = TRUE)
rac_knots <- list(month = c(0.5, 12.5))
rodent_m <- bam(abundance ~
  s(plot_f, bs = "re") +
  t2(
    rank,
    month,
    year,
    by = plot_f,
    bs = c("cr", "cc", "cr"),
    k = c(5, 12, 20),
    full = TRUE
  ),
data = rodent_rac,
family = poisson(),
method = "fREML",
discrete = TRUE,
nthreads = 3,
control = ctrl,
knots = rac_knots
)

ds <- data_slice(
  rodent_m,
  rank = evenly(rank, n = 50),
  month = evenly(month, by = 1),
  year = evenly(year, n = 50),
  plot_f = evenly(plot_f)
)

fv_m1 <- fitted_values(rodent_m, data = ds)

fv_m1 |>
  filter(plot_f == "4") |>
  ggplot(aes(x = rank, y = .fitted, colour = year, group = year)) +
  geom_line() +
  facet_wrap(~ month)

# Plot 4
rodent_m_p4 <- bam(abundance ~
  # t2(
  #   rank,
  #   month,
  #   year,
  #   bs = c("cr", "cc", "cr"),
  #   k = c(5, 12, 20),
  #   full = TRUE
  # ),
  s(rank, bs = "cr", k = 6) +
  s(rank, census_f, bs = "fs", k = 6),
data = rodent_rac |> filter(plot_f == "17"),
family = nb(),
method = "fREML",
discrete = TRUE,
nthreads = 3,
control = ctrl,
knots = rac_knots
)

rodent_rac |>
  filter(plot_f == "4") |>
  ggplot(aes(x = rank, y = abundance,
    colour = year, group = censusdate)) +
  geom_line(alpha = 0.3) +
  scale_colour_viridis_c(option = "magma")

ds <- data_slice(
  rodent_m_p4,
  rank = evenly(rank, n = 25),
  census_f = evenly(census_f)
)

fv_p4 <- fitted_values(rodent_m_p4, data = ds)

fv_p4 |>
  ggplot(aes(x = rank, y = .fitted, colour = year, group = year)) +
  geom_line() +
  facet_wrap(~ month)