# assign("pkg.globals", new.env())

.onLoad <- function(libname, pkgname)
{
  # >>>
  # assign("pkg.globals", new.env())
  assign('density_control_rds',
         system.file("extdata",
                     "healthycontrol.fragmentprofile.RDS",
                     package = "cfdnakit"), envir = .GlobalEnv)
  # >>> Set default control for distribution plot
  # pkg.globals$density_control_rds =
    # system.file("extdata",
    #             "healthycontrol.fragmentprofile.RDS",
    #             package = "cfdnakit")
  # <<<
  # <<<
}
