library(here)

path_default_nix <- here() 

library(rix)

rix(date = "2025-09-09",
    r_pkgs = c("languageserver", "rix", "here", "pacman", "data.table", "dplyr", "ggplot2", "fixest", "texreg", "kableExtra", "broom", "glue", "modelsummary"),
    system_pkgs = c("quarto", "git"),
    git_pkgs = NULL,
    ide = "none",
    project_path = path_default_nix,
    overwrite = TRUE,
    print = TRUE)
