library(tidyverse)
list.files("outputData/SSOM/marginal_effect")

marginal_dir <- "outputData/SSOM/marginal_effect"

files <- list.files(
  marginal_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

marginal_all <- map_dfr(files, function(f) {
  
  dat <- read_csv(f, show_col_types = FALSE)
  fname <- basename(f)
  
  tibble(
    covariate = str_extract(fname, "(?<=SSOM_effect_)[A-Za-z]+"),
    Year      = str_extract(fname, "\\d{4}")
  ) %>%
    bind_cols(dat)
})


marginal_all <- marginal_all %>%
  mutate(
    Year = factor(Year, levels = sort(unique(Year))),
    covariate = factor(covariate,
                       levels = c("AG", "ndvi", "Rf", "Tmax"))
  )

library(tidyverse)
library(patchwork)

plot_marginal <- function(df, yr, cov) {
  
  ggplot(
    df %>% filter(Year == yr, covariate == cov),
    aes(x = x, y = fit)
  ) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                alpha = 0.25, fill = "grey60") +
    geom_line(size = 1) +
    labs(
      title = paste( yr),
      #title = paste(cov, "-", yr),
      
      x = cov,
      y = "Occupancy (ψ)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text  = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12)
    )
}

years <- sort(unique(marginal_all$Year))
covs  <- levels(marginal_all$covariate)

all_plots <- list()

k <- 1
for (yr in years) {
  for (cv in covs) {
    
    if (nrow(filter(marginal_all, Year == yr, covariate == cv)) == 0)
      next
    
    all_plots[[k]] <- plot_marginal(marginal_all, yr, cv)
    k <- k + 1
  }
}


plots_per_page <- 12

plot_pages <- split(
  all_plots,
  ceiling(seq_along(all_plots) / plots_per_page)
)


for (i in seq_along(plot_pages)) {
  print(
    wrap_plots(plot_pages[[i]], ncol = 4, nrow = 3)
  )
}


for (i in seq_along(plot_pages)) {
  
  ggsave(
    filename = paste0(
      "outputData/SSOM/figures/marginal_effects_page_", i, ".png"
    ),
    plot   = wrap_plots(plot_pages[[i]], ncol = 4, nrow = 3),
    width  = 9,
    height = 8.5,
    units  = "in",
    dpi    = 300
  )
}


pdf(
  file   = "outputData/SSOM/figures/marginal_effects_A4_portrait.pdf",
  width  = 8.27,
  height = 11.69,
  onefile = TRUE
)

for (i in seq_along(plot_pages)) {
  print(
    wrap_plots(plot_pages[[i]], ncol = 4, nrow = 3)
  )
}

dev.off()
