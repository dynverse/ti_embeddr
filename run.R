#!/usr/local/bin/Rscript

task <- dyncli::main()

# load libraries
library(dyncli, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

library(parallel, warn.conflicts = FALSE)
library(BiocGenerics, warn.conflicts = FALSE)
library(Biobase, warn.conflicts = FALSE)

library(scaterlegacy, warn.conflicts = FALSE)
library(embeddr, warn.conflicts = FALSE)

#####################################
###           LOAD DATA           ###
#####################################
expression <- task$expression
params <- task$params
priors <- task$priors

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################
# calculate nn param
nn <- max(round(log(nrow(expression)) * params$nn_pct), 9)

# load data in scaterlegacy
sce <- scaterlegacy::newSCESet(exprsData = t(expression))

# run embeddr
sce <- embeddr::embeddr(
  sce,
  kernel = params$kernel,
  metric = params$metric,
  nn = nn,
  eps = params$eps,
  t = params$t,
  symmetrize = params$symmetrize,
  measure_type = params$measure_type,
  p = params$ndim
)

# fit pseudotime
sce <- embeddr::fit_pseudotime(
  sce,
  thresh = params$thresh,
  maxit = params$maxit,
  stretch = params$stretch,
  smoother = params$smoother
)

# TIMING: done with trajectory inference
timings$method_aftermethod <- Sys.time()

#####################################
###     SAVE OUTPUT TRAJECTORY    ###
#####################################

# construct milestone network
pseudotime <-
  as(sce@phenoData, "data.frame")$pseudotime %>%
  setNames(rownames(expression))

# creating extra output for visualisation purposes
dimred <- sce@reducedDimension

dimred_segment_points <-
  as(sce@phenoData, "data.frame") %>%
  dplyr::arrange(pseudotime) %>%
  select(starts_with("trajectory_")) %>%
  as.matrix()

dimred_segment_progressions <-
  tibble(
    from = "milestone_begin",
    to = "milestone_end",
    percentage = pseudotime %>% sort
  )

output <-
  wrap_data(
    cell_ids = rownames(expression)
  ) %>%
  add_linear_trajectory(
    pseudotime = pseudotime,
    do_scale_minmax = FALSE
  ) %>%
  add_dimred(
    dimred = dimred,
    dimred_milestones = dimred_segment_points[c(1, nrow(dimred_segment_points)), , drop = FALSE],
    dimred_segment_points = dimred_segment_points,
    dimred_segment_progressions = dimred_segment_progressions
  ) %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
