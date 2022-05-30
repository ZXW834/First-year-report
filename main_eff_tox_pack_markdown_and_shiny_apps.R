source("calcp_shiny_app.R")
shinyApp(ui,server)

if (pval>0){
  rmarkdown::render("D:/University of Southampton/PhD Supervisor Dave and Dankmar/R/Rmarkdown/rstan/rstan/Eff-Tox_design_stan.Rmd")
}

