library(here)

path_default_nix <- here() 

library(rix)

rix(date = "2025-11-10",
    r_pkgs = c("languageserver", "rix", "here", "pacman", "data.table", "dplyr", "ggplot2", "fixest", "texreg", "kableExtra", "broom", "glue", "modelsummary"),
    system_pkgs = c("quarto", "git"),
    git_pkgs = NULL,
    jl_conf = list(
		   jl_version = "1.11",
		   jl_pkgs = c("DataFrames", "CSV", "LinearAlgebra", "Statistics")
		   ),
    ide = "none",
    project_path = path_default_nix,
    overwrite = TRUE,
    print = TRUE)
