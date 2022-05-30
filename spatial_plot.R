library(tidyverse)

# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
Time.step <- 10

### IID
df <- data.frame()

for (Nparticles in c(100)) {
  df_dac <- read.csv(paste0("data/spatial/stats_dac_spatial_tau", tau, "d", d, "N", Nparticles),
                      row.names = 1)
  mse_dac <- df_dac[1:d, 1:d]
  colnames(mse_dac) <- paste0("c",seq(1,d))
  rownames(mse_dac) <- paste0("r",seq(1,d))
  tmp <- mse_dac %>%
    as.data.frame() %>%
    rownames_to_column("r_id") %>%
    pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
    mutate(col = fct_relevel(col, colnames(mse_dac)))
  tmp$N <- Nparticles
  tmp$runtime <- df_dac$runtime[1]
  df <- rbind(df, tmp)
}

ggplot(data = df, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = mse)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~N, ncol = 3, nrow = 1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank(), panel.spacing = unit(1, "lines"),
        legend.title = element_blank(), legend.text=element_text(size=25),
        text = element_text(size=20), aspect.ratio = 1)
# ggsave("spatial_mse.pdf", width = 20, height = 15, dpi = 300)
