library(tidyverse)

# dimension
d <- 8
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 100

Nparticles <- 1000
M <- 100

# nsmc
df_nsmc <- read.csv(paste0("data/nl/stats_nsmc_nl_iid_d", d, "N", Nparticles),
                    row.names = 1)
mse_nsmc <- df_nsmc[1:d, 1:d]
colnames(mse_nsmc) <- paste0("c",seq(1,d))
rownames(mse_nsmc) <- paste0("r",seq(1,d))
var_nsmc <- df_nsmc[(d+1):(2*d), 1:d]
colnames(var_nsmc) <- paste0("c",seq(1,d))
rownames(var_nsmc) <- paste0("r",seq(1,d))
df1 <- mse_nsmc %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
  mutate(col = fct_relevel(col, colnames(mse_nsmc)))
df1$algo <- "nsmc"
df1var <- var_nsmc %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "var") %>%
  mutate(col = fct_relevel(col, colnames(var_nsmc)))
df1var$algo <- "nsmc"

# stpf
df_stpf <- read.csv(paste0("data/nl/stats_stpf_nl_iid_d", d, "N", Nparticles),
                    row.names = 1)
mse_stpf <- df_stpf[1:d, 1:d]
colnames(mse_stpf) <- paste0("c",seq(1,d))
rownames(mse_stpf) <- paste0("r",seq(1,d))
df2 <- mse_stpf %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
  mutate(col = fct_relevel(col, colnames(mse_stpf)))
df2$algo <- "stpf"
var_stpf <- df_stpf[(d+1):(2*d), 1:d]
colnames(var_stpf) <- paste0("c",seq(1,d))
rownames(var_stpf) <- paste0("r",seq(1,d))
df2var <- var_stpf %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "var") %>%
  mutate(col = fct_relevel(col, colnames(var_stpf)))
df2var$algo <- "stpf"

# dac
df_dac <- read.csv(paste0("data/nl/stats_dac_nl_iid_d", d, "N", Nparticles),
                    row.names = 1)
mse_dac <- df_dac[1:d, 1:d]
colnames(mse_dac) <- paste0("c",seq(1,d))
rownames(mse_dac) <- paste0("r",seq(1,d))
df3 <- mse_dac %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
  mutate(col = fct_relevel(col, colnames(mse_dac)))
df3$algo <- "dac"
var_dac <- df_dac[(d+1):(2*d), 1:d]
colnames(var_dac) <- paste0("c",seq(1,d))
rownames(var_dac) <- paste0("r",seq(1,d))
df3var <- var_dac %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "var") %>%
  mutate(col = fct_relevel(col, colnames(var_dac)))
df3var$algo <- "dac"

# MSE
df <- rbind(df1, df2, df3)

ggplot(data = df, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = mse)) +
geom_raster() +
scale_fill_viridis_c() +
facet_wrap(~algo, ncol = 2, nrow = 2) +
theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15), strip.text.x = element_text(size = 20))
# ggsave("nl_d32N100_mse.pdf", width = 12, height = 12, dpi = 300)

# variance
df_var <- rbind(df1var, df2var, df3var)
ggplot(data = df_var, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = var)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~algo, ncol = 2, nrow = 2) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15), strip.text.x = element_text(size = 20))
# ggsave("nl_d32N100_var.pdf", width = 12, height = 12, dpi = 300)
