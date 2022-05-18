library(tidyverse)

# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 10

### IID
df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

for (Nparticles in c(100, 1000)) {
  # nsmc
  df_nsmc <- read.csv(paste0("data/nl/stats_nsmc_nl_iid_d", d, "N", Nparticles),
                      row.names = 1)
  mse_nsmc <- df_nsmc[1:d, 1:d]
  colnames(mse_nsmc) <- paste0("c",seq(1,d))
  rownames(mse_nsmc) <- paste0("r",seq(1,d))
  tmp <- mse_nsmc %>%
    as.data.frame() %>%
    rownames_to_column("r_id") %>%
    pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
    mutate(col = fct_relevel(col, colnames(mse_nsmc)))
  tmp$N <- Nparticles
  tmp$runtime <- df_nsmc$runtime[1]
  print(paste(tmp$runtime[1]))
  df1 <- rbind(df1, tmp)
  # stpf
  df_stpf <- read.csv(paste0("data/nl/stats_stpf_nl_iid_d", d, "N", Nparticles),
                      row.names = 1)
  mse_stpf <- df_stpf[1:d, 1:d]
  colnames(mse_stpf) <- paste0("c",seq(1,d))
  rownames(mse_stpf) <- paste0("r",seq(1,d))
  tmp <- mse_stpf %>%
    as.data.frame() %>%
    rownames_to_column("r_id") %>%
    pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
    mutate(col = fct_relevel(col, colnames(mse_stpf)))
  tmp$N <- Nparticles
  tmp$runtime <- df_stpf$runtime[1]
  print(paste(tmp$runtime[1]))
  df2 <- rbind(df2, tmp)
  # dac
  df_dac <- read.csv(paste0("data/nl/stats_dac_nl_iid_d", d, "N", Nparticles),
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
  print(paste(tmp$runtime[1]))
  df3 <- rbind(df3, tmp)
}
df1$algo <- "nsmc"
df2$algo <- "stpf"
df3$algo <- "dac"

# MSE
df <- rbind(df1, df2, df3)
df$type <- "iid"
ggplot(data = df, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = mse)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~interaction(algo, N), ncol = 3, nrow = 2) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank(), panel.spacing = unit(1, "lines"),
        legend.title = element_blank(), legend.text=element_text(size=25),
        text = element_text(size=20), aspect.ratio = 1)
# ggsave("nl_mse.pdf", width = 20, height = 15, dpi = 300)

df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()
### NON IID
for (Nparticles in c(100)) {
  # nsmc
  df_nsmc <- read.csv(paste0("data/nl/stats_nsmc_nl_cov_d", d, "N", Nparticles),
                      row.names = 1)
  mse_nsmc <- df_nsmc[1:d, 1:d]
  colnames(mse_nsmc) <- paste0("c",seq(1,d))
  rownames(mse_nsmc) <- paste0("r",seq(1,d))
  tmp <- mse_nsmc %>%
    as.data.frame() %>%
    rownames_to_column("r_id") %>%
    pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
    mutate(col = fct_relevel(col, colnames(mse_nsmc)))
  tmp$N <- Nparticles
  tmp$runtime <- df_nsmc$runtime[1]
  df1 <- rbind(df1, tmp)
  # stpf
  df_stpf <- read.csv(paste0("data/nl/stats_stpf_nl_cov_d", d, "N", Nparticles),
                      row.names = 1)
  mse_stpf <- df_stpf[1:d, 1:d]
  colnames(mse_stpf) <- paste0("c",seq(1,d))
  rownames(mse_stpf) <- paste0("r",seq(1,d))
  tmp <- mse_stpf %>%
    as.data.frame() %>%
    rownames_to_column("r_id") %>%
    pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
    mutate(col = fct_relevel(col, colnames(mse_stpf)))
  tmp$N <- Nparticles
  tmp$runtime <- df_stpf$runtime[1]
  df2 <- rbind(df2, tmp)
  # dac
  df_dac <- read.csv(paste0("data/nl/stats_dac_nl_cov_d", d, "N", Nparticles),
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
  df3 <- rbind(df3, tmp)
}
df1$algo <- "nsmc"
df2$algo <- "stpf"
df3$algo <- "dac"

# MSE
df_cov <- rbind(df1, df2, df3)
df_cov$type <- "non-iid"
ggplot(data = df_cov, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = mse)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~interaction(algo, N), ncol = 3, nrow = 2) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank(), panel.spacing = unit(1, "lines"),
        legend.title = element_blank(), legend.text=element_text(size=25),
        text = element_text(size=20), aspect.ratio = 1)
# ggsave("nl_mse_cov.pdf", width = 20, height = 15, dpi = 300)
df_collective <- aggregate(cbind(mse, runtime) ~ algo + N + type, data = rbind(df, df_cov), FUN = "mean")
ggplot(data = df_collective, aes(x=runtime, y=mse, group = interaction(algo, type), color = algo))+
  geom_line(aes(linetype = as.factor(type)))+
  geom_point(aes(shape = as.factor(N)))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
