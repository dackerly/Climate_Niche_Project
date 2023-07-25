options(scipen=999)

# load in processed FIA data from Rosenblad et al. (2023) PNAS
# niche means were calculated from western continental US, not just CA
data <- readRDS("./data/kr_fia_data_products/subplot_data.RDS")
data <- subset(data, timepoint=="after" & STATECD==6 & !is.na(STATECD))

library(ggplot2)
opt_meansimple <- ggplot(data, aes(x=tmp_nicheopt_after/10, y=tmp_nichemeansimple_after/10))+
  geom_point()+
  xlab("Modeled Niche Optimum (°C)")+
  ylab("Basal Area-Weighted Mean of Presences (°C)")+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method=lm)+
  theme_bw()
opt_meansimple
ggsave("./results/kr_fia_results/fia_niches_opt_meansimple.png", height=5, width=5)

opt_mean <- ggplot(data, aes(x=tmp_nicheopt_after/10, y=tmp_nichemean_after/10))+
  geom_point()+
  xlab("Modeled Niche Optimum (°C)")+
  ylab("Modeled Niche Mean (°C)")+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method=lm)+
  theme_bw()
opt_mean
ggsave("./results/kr_fia_results/fia_niches_opt_mean.png", height=5, width=5)

opt_mse_1to1 <- mean((data$tmp_nicheopt_after-data$tmp)^2)
mean_mse_1to1 <- mean((data$tmp_nichemean_after-data$tmp)^2)
meansimple_mse_1to1 <- mean((data$tmp_nichemeansimple_after-data$tmp)^2)

opt_clim_mod <- lm(tmp_nicheopt_after ~ tmp, data=data)
mean_clim_mod <- lm(tmp_nichemean_after ~ tmp, data=data)
meansimple_clim_mod <- lm(tmp_nichemeansimple_after ~ tmp, data=data)

opt_r2 <- summary(opt_clim_mod)$adj.r.squared
mean_r2 <- summary(mean_clim_mod)$adj.r.squared
meansimple_r2 <- summary(meansimple_clim_mod)$adj.r.squared

opt_es <- opt_clim_mod$coefficients[2]
mean_es <- mean_clim_mod$coefficients[2]
meansimple_es <- meansimple_clim_mod$coefficients[2]

opt_clim <- ggplot(data, aes(x=tmp/10, y=tmp_nicheopt_after/10))+
  geom_point(pch=".")+
  xlab("Mean Annual Temperature (°C)")+
  ylab("Modeled Niche Optimum (°C)")+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method=lm)+
  theme_bw()
opt_clim
ggsave("./results/kr_fia_results/fia_niches_opt_clim.png", height=5, width=5)

mean_clim <- ggplot(data, aes(x=tmp/10, y=tmp_nichemean_after/10))+
  geom_point(pch=".")+
  xlab("Mean Annual Temperature (°C)")+
  ylab("Modeled Niche Mean (°C)")+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method=lm)+
  theme_bw()
mean_clim
ggsave("./results/kr_fia_results/fia_niches_mean_clim.png", height=5, width=5)

meansimple_clim <- ggplot(data, aes(x=tmp/10, y=tmp_nichemeansimple_after/10))+
  geom_point(pch=".")+
  xlab("Mean Annual Temperature (°C)")+
  ylab("Basal Area-Weighted Mean of Presences (°C)")+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method=lm)+
  theme_bw()
meansimple_clim
ggsave("./results/kr_fia_results/fia_niches_meansimple_clim.png", height=5, width=5)

summ <- data.frame(metric=c("opt", "mean", "meansimple"),
                   mse1to1=c(opt_mse_1to1, mean_mse_1to1, meansimple_mse_1to1),
                   r2=c(opt_r2, mean_r2, meansimple_r2),
                   es=c(opt_es, mean_es, meansimple_es))

ggplot(summ, aes(x=metric, y=es, color=r2, size=mse1to1))+
  ggtitle("MAT")+
  geom_hline(yintercept=1, lty="dotted")+
  geom_point()+
  xlab("Community Weighted Niche Metric")+
  ylab("Effect Size of Macroclimatic Predictor")+
  scale_x_discrete(labels=c("mean"="Modeled\nMean", "opt"="Modeled\nOptimum", "meansimple"="Basal Area-Weighted\nMean of Presences"))+
  scale_color_viridis_c(name="R²")+
  scale_size_continuous(name="1:1\nMSE",
                        range=c(
                          min(summ$mse1to1)/100,
                          max(summ$mse1to1)/100))+
  theme_bw()
ggsave("./results/kr_fia_results/fia_niches_and_plots_summary_tmp.png", height=5, width=5)
