library(tidyverse)

# dimension
d <- 32
# parameters
sigmaX <- 1
nu <- 10
tau <- 1/4
delta <- 1
Time.step <- 10

Nparticles <- 100
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

# dac wit covariance
df_dac_cov <- read.csv(paste0("data/nl/stats_dac_nl_cov_d", d, "N", Nparticles),
                   row.names = 1)
mse_dac_cov <- df_dac_cov[1:d, 1:d]
colnames(mse_dac_cov) <- paste0("c",seq(1,d))
rownames(mse_dac_cov) <- paste0("r",seq(1,d))
df4 <- mse_dac_cov %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "mse") %>%
  mutate(col = fct_relevel(col, colnames(mse_dac_cov)))
df4$algo <- "dac-cov"
var_dac_cov <- df_dac_cov[(d+1):(2*d), 1:d]
colnames(var_dac_cov) <- paste0("c",seq(1,d))
rownames(var_dac_cov) <- paste0("r",seq(1,d))
df4var <- var_dac_cov %>%
  as.data.frame() %>%
  rownames_to_column("r_id") %>%
  pivot_longer(-c(r_id), names_to = "col", values_to = "var") %>%
  mutate(col = fct_relevel(col, colnames(var_dac_cov)))
df4var$algo <- "dac-cov"

# MSE
df <- rbind(df1, df2, df3, df4)

ggplot(data = df, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = mse)) +
geom_raster() +
scale_fill_viridis_c() +
facet_wrap(~algo, ncol = 1, nrow = 4) +
theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15), strip.text.x = element_text(size = 20))
# ggsave("nl_d32N100_mse.pdf", width = 5, height = 15, dpi = 300)

# variance
df_var <- rbind(df1var, df2var, df3var, df4var)
ggplot(data = df_var, aes(x = col , y = factor(r_id, level = paste0("r", d:1)), fill = var)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~algo, ncol = 1, nrow = 4) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_blank(), legend.text=element_text(size=20),
        text = element_text(size=15), strip.text.x = element_text(size = 20))
# ggsave("nl_d32N100_var.pdf", width = 5, height = 15, dpi = 300)

df_nsmc2 <- read.csv(paste0("data/nl/stats_nsmc_nl_iid_d", d, "N", 1000),
                    row.names = 1)

df_collective <- data.frame(c(c(df_stpf$runtime[1], df_dac$runtime[1], df_nsmc$runtime[1])/10, df_dac_cov$runtime[1]),
                            c(mean(as.matrix(mse_stpf)), mean(as.matrix(mse_dac)), mean(as.matrix(mse_nsmc)), mean(as.matrix(mse_dac_cov))),
                            c("stpf", "dac", "nsmc", "dac-cov"))
colnames(df_collective) <- c("runtime", "mse", "algo")
df_collective$N <- "10^2"
df_collective2 <- data.frame(c(df_nsmc2$runtime[1]),
                            c(mean(as.matrix(df_nsmc2[1:d, 1:d]))),
                            c("nsmc"))
colnames(df_collective2) <- c("runtime", "mse", "algo")
df_collective2$N <- "10^3"
df_collective <- rbind(df_collective, df_collective2)
ggplot(data = df_collective, aes(x=runtime, y=mse, group = algo, color = algo))+
  geom_line()+
  geom_point(aes(shape = N))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
