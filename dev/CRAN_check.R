# Build and check ernm for CRAN submission.
#
# Usage:
#   Rscript dev/CRAN_check.R check        # document, build, R CMD check --as-cran
#   Rscript dev/CRAN_check.R build        # document and build the source tarball only
#   Rscript dev/CRAN_check.R diagnostics  # spelling, URL, and roxygen tag checks
#   Rscript dev/CRAN_check.R revdep       # optional reverse dependency check
#   Rscript dev/CRAN_check.R all          # check plus diagnostics
#
# By default, CRAN check fails on ERROR, WARNING, or NOTE. Set
# ERNM_CRAN_STRICT=false to allow NOTE output while still failing on ERROR/WARNING.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "check" else args[[1]]

usage <- function() {
  cat(
    "Usage: Rscript dev/CRAN_check.R [check|build|diagnostics|revdep|all|help]\n",
    "\n",
    "Modes:\n",
    "  check        document, build, and run R CMD check --as-cran\n",
    "  build        document and build the CRAN source tarball only\n",
    "  diagnostics  run optional spelling, URL, and roxygen tag checks\n",
    "  revdep       run optional reverse dependency checks\n",
    "  all          run check and diagnostics\n",
    "\n",
    "Environment:\n",
    "  ERNM_CRAN_STRICT=false  allow CRAN check NOTE output\n",
    sep = ""
  )
}

if (mode %in% c("help", "--help", "-h")) {
  usage()
  quit(status = 0)
}

valid_modes <- c("check", "build", "diagnostics", "revdep", "all")
if (!mode %in% valid_modes) {
  usage()
  stop("Unknown CRAN check mode: ", mode, call. = FALSE)
}

script_dir <- function() {
  file_arg <- "--file="
  command_args <- commandArgs(trailingOnly = FALSE)
  script_arg <- command_args[startsWith(command_args, file_arg)][1]
  if (!is.na(script_arg)) {
    return(dirname(normalizePath(sub(file_arg, "", script_arg), mustWork = TRUE)))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE)))
  }

  getwd()
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Install required package before running this script: install.packages(\"",
      package,
      "\")",
      call. = FALSE
    )
  }
}

require_optional_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("Skipping optional check; package is not installed: ", package)
    return(FALSE)
  }
  TRUE
}

description_packages <- function(description, fields) {
  values <- unname(description[intersect(fields, names(description))])
  if (length(values) == 0) {
    return(character())
  }

  packages <- unlist(strsplit(paste(values, collapse = ","), ",", fixed = TRUE))
  packages <- sub("\\s*\\(.*\\)", "", packages)
  packages <- trimws(packages)
  setdiff(packages[nzchar(packages)], "R")
}

check_cran_dependencies <- function(description) {
  packages <- description_packages(description, c("Depends", "Imports", "Suggests", "LinkingTo"))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Install package dependencies before running CRAN check: install.packages(c(",
      paste(sprintf("\"%s\"", missing), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
}

load_source_with_data <- function(path) {
  env <- roxygen2::load_source(path)
  data_dir <- file.path(path, "data")
  if (dir.exists(data_dir)) {
    data_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
    lapply(data_files, load, envir = env)
  }
  env
}

document_package <- function(package_path) {
  require_package("roxygen2")

  message("Documenting package without loading compiled Rcpp modules...")
  roxygen2::roxygenise(package_path, load_code = load_source_with_data)
}

build_cran_tarball <- function(package_path, build_dir) {
  require_package("devtools")

  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  document_package(package_path)

  message("Building CRAN source tarball...")
  devtools::build(
    pkg = package_path,
    path = build_dir,
    vignettes = TRUE,
    manual = TRUE,
    quiet = FALSE
  )
}

check_cran_tarball <- function(tarball, check_dir) {
  require_package("rcmdcheck")

  dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

  old_repos <- getOption("repos")
  on.exit(options(repos = old_repos), add = TRUE)
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  message("Running R CMD check --as-cran...")
  result <- rcmdcheck::rcmdcheck(
    path = tarball,
    args = c("--as-cran"),
    check_dir = check_dir,
    error_on = "never"
  )

  print(result)

  strict <- !tolower(Sys.getenv("ERNM_CRAN_STRICT", "true")) %in% c("false", "0", "no")
  has_errors <- length(result$errors) > 0
  has_warnings <- length(result$warnings) > 0
  has_notes <- length(result$notes) > 0

  if (has_errors || has_warnings || (strict && has_notes)) {
    stop(
      "CRAN check did not pass cleanly: ",
      length(result$errors), " error(s), ",
      length(result$warnings), " warning(s), ",
      length(result$notes), " note(s).",
      call. = FALSE
    )
  }

  invisible(result)
}

run_diagnostics <- function(package_path) {
  message("Running optional CRAN diagnostics...")

  if (require_optional_package("checkhelper")) {
    checkhelper::find_missing_tags(path = package_path)
  }

  if (require_optional_package("spelling")) {
    spelling::spell_check_package(package_path)
  }

  if (require_optional_package("urlchecker")) {
    urlchecker::url_check(package_path)
  }

  invisible(TRUE)
}

run_revdep <- function(package_path, revdep_dir) {
  require_package("revdepcheck")

  dir.create(revdep_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running reverse dependency checks...")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(package_path)

  revdepcheck::revdep_check(
    pkg = package_path,
    revdep_dir = revdep_dir,
    num_workers = max(1, parallel::detectCores() - 1)
  )
  revdepcheck::revdep_summary(revdep = revdep_dir)
  revdepcheck::revdep_report(revdep = revdep_dir)
}

package_path <- normalizePath(file.path(script_dir(), ".."), mustWork = TRUE)
project_dir <- dirname(package_path)
description <- read.dcf(file.path(package_path, "DESCRIPTION"))[1, ]
package_name <- unname(description[["Package"]])
package_version <- unname(description[["Version"]])
build_dir <- file.path(project_dir, "cran-build")
check_dir <- file.path(project_dir, "cran-check")
revdep_dir <- file.path(project_dir, "revdep")

message("Package: ", package_name, " ", package_version)
message("Package path: ", package_path)

if (mode %in% c("build", "check", "all")) {
  if (mode %in% c("check", "all")) {
    check_cran_dependencies(description)
  }

  tarball <- build_cran_tarball(package_path, build_dir)
  message("Built source tarball: ", tarball)
}

if (mode %in% c("check", "all")) {
  check_cran_tarball(tarball, check_dir)
}

if (mode %in% c("diagnostics", "all")) {
  run_diagnostics(package_path)
}

if (mode == "revdep") {
  run_revdep(package_path, revdep_dir)
}

message("Done.")
