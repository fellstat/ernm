# Prepare for CRAN ----
# adatped from : https://github.com/ThinkR-open/prepare-for-cran

# no notes no warnings no errors
devtools::check()

# Update dependencies in DESCRIPTION
# install.packages('attachment', repos = 'https://thinkr-open.r-universe.dev')
# attachment::att_amend_desc()

# Check package coverage
# covr::package_coverage()
# covr::report()

# Run tests
devtools::test()

# Check package as CRAN using the correct CRAN repo
withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org/")),
                    {callr::default_repos()
                      rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran")) })

# Check content
checkhelper::find_missing_tags()

# Check spelling - No typo
# usethis::use_spell_check()
spelling::spell_check_package()

# Check URL are correct
# install.packages('urlchecker', repos = 'https://r-lib.r-universe.dev')
urlchecker::url_check()
urlchecker::url_update()

# check on other distributions
# _rhub v2
rhub::rhub_setup() # Commit, push, merge
rhub::rhub_doctor()
rhub::rhub_platforms()
plats <- rhub::rhub_platforms()
no_run <- grepl('unstable',plats$r_version)
for(plat in plats[[1]][!no_run]){
  rhub::rhub_check(platform = plat)
}

# Check reverse dependencies
#remotes::install_github("r-lib/revdepcheck")
usethis::use_git_ignore("revdep/")
usethis::use_build_ignore("revdep/")

devtools::revdep()
library(revdepcheck)
# In another session because Rstudio interactive change your config:
id <- rstudioapi::terminalExecute("Rscript -e 'revdepcheck::revdep_check(num_workers = 4)'")
rstudioapi::terminalKill(id)
# if [Exit Code] is not 0, there is a problem !
# to see the problem: execute the command in a new terminal manually.

# See outputs now available in revdep/
revdep_details(revdep = "pkg")
revdep_summary()                 # table of results by package
revdep_report()
# Clean up when on CRAN
revdep_reset()

# Verify you're ready for release, and release
if(FALSE){
    # make comments for CRAN
    usethis::use_cran_comments(open = rlang::is_interactive())
    # Upgrade version number
    usethis::use_version(which = c("patch", "minor", "major", "dev")[1])
    devtools::release()
}