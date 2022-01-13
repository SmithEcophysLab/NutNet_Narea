# Script using least-cost to explore Narea response at nutnet
## generally follows the analysis structure of Dong et al. (2017) Biogeosciences
## but includes the addition of soil N as a predictor

#### load packages ####
library(tidyverse)
library(lme4)
library(car)
library(r2glmm)
library(treemapify)
library(emmeans)
library(relaimpo)
library(patchwork)
library(multcomp)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)

#### load functions ####
### functions to calculate vcmax and jmax 
source('optimal_vcmax_R/calc_optimal_vcmax.R')
sourceDirectory('optimal_vcmax_R/functions', modifiedOnly = FALSE)

### function to calculate relative importance for mixed models 
### from https://gist.github.com/BERENZ/e9b581a4b7160357934e
calc.relip.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  importances <- calc.relimp(as.matrix(cov_XY), rela = F, type = type)
  return(importances)
}

calc.relip.boot.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  bootresults <- boot.relimp(as.matrix(cov_XY), b=1000, rela = F, type = type)
  importances <- booteval.relimp(bootresults, norank=T)
  return(importances)
}


#### hypothesis figure ####
### leaf N and chi by N supply with different differences in N demand

n_supply_trend <- calc_optimal_vcmax(beta = seq(100, 300, 10))
n_supply_trend$photo_n_nobetachange <- n_supply_trend$nphoto[10]

hypothesis_data <- data.frame(cbind(c(n_supply_trend$beta, n_supply_trend$beta), 
                        c(n_supply_trend$nphoto, n_supply_trend$photo_n_nobetachange),
                        c(rep('no_change', 21), rep('change', 21))))
colnames(hypothesis_data) <- c('beta', 'Narea', 'demand')
hypothesis_data$beta <- as.numeric(as.character(hypothesis_data$beta))
hypothesis_data$Narea <- as.numeric(as.character(hypothesis_data$Narea))

(hypothesis_plot <- ggplot(data = hypothesis_data, 
                         aes(x = (1/beta), y = Narea, col = demand)) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 25),
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(size = 4, aes(linetype = demand, color = demand)) +
  scale_linetype_manual(values = c(1, 2), guide = "none") +
  scale_colour_manual(values = c('black', 'grey'), 
                      labels = c('∆AGB = high', '∆AGB = low')) +
  guides(color = guide_legend(title = NULL)) +
  ylab(expression(italic('N')['area'])) +
  xlab(expression(italic('N')['availability'])))

#### map of sites ####
site_data <- read.csv("../Data/site_metadata.csv")
sites <- site_data %>% 
  subset(site_code == "bldr.us" | site_code == "bogong.au" | 
           site_code == "burrawan.au" | site_code == "cbgb.us" | site_code == "comp.pt" | 
           site_code == "cowi.ca" | site_code == "elliot.us" | 
           site_code == "gilb.za" | site_code == "hopl.us" | site_code == "lancaster.uk" | 
           site_code == "mcla.us" | site_code == "mtca.au" | 
           site_code == "sage.us" | site_code == "sgs.us" | 
           site_code == "sier.us" | site_code == "smith.us" | 
           site_code == "summ.za" | site_code == "unc.us" | site_code == "valm.ch")

world_map <- map_data("world")
world_map_subset <- subset(world_map, region != "Antarctica")
(map_plot <- ggplot() + geom_polygon(data = world_map_subset, 
                                     aes(x = long, y = lat, group = group),
                                     fill = 'white', color = 'gray35') + 
    coord_fixed(1.3) +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = 'white', colour = NA)) +
    geom_point(data = sites, aes(x = longitude, y = latitude), 
               shape = 16, color = "red", size = 1) +
    geom_point(data = sites, aes(x = longitude, y = latitude), 
               shape = 1, color = "red", fill = NA, size = 3, stroke = 1))

#### drivers of Narea and their relative importance ####
### load data
leaf <- read.csv('../Data/leaf_traits.csv')

## turn treatment numbers into factors
leaf$Ntrt_fac <- as.factor(leaf$Ntrt)
leaf$Ptrt_fac <- as.factor(leaf$Ptrt)
leaf$Ktrt_fac <- as.factor(leaf$Ktrt)
leaf$block_fac <- as.factor(leaf$block)

## scale temperature data
leaf$tmp_scaled <- leaf$tmp - 25

## assign plant functional groups
leaf$pft <- leaf$functional_group
leaf$pft[leaf$pft == 'NULL'] <- NA
leaf$grass <- 'no'
leaf$grass[leaf$pft == 'GRASS'] <- 'yes'
leaf$fence <- 'no'
leaf$fence[leaf$trt == 'Fence' | leaf$trt == 'NPK+Fence'] <- 'yes'

## create subset of data where chi is greater than 0 and less than 1
leaf_chi_subset = subset(leaf, chi > 0 & chi < 1) # lose 659 points

### linear mixed effects model
leaf_chi_subset$logpar <- log(leaf_chi_subset$par)
leaf_chi_subset$logpar_per_leaf_area <- log(leaf_chi_subset$par_per_leaf_area)
leaf_chi_subset$logvpd <- log(leaf_chi_subset$vpd)
leaf_chi_subset$loglma <- log(leaf_chi_subset$lma)
leaf$logpar <- log(leaf$par)
leaf$logpar_per_leaf_area <- log(leaf$par_per_leaf_area)
leaf$logvpd <- log(leaf$vpd)
leaf$loglma <- log(leaf$lma)
leafNarea_lmer <- lmer(log(narea) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + 
                         logpar_per_leaf_area + loglma + chi + Nfix + photosynthetic_pathway +
                        (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                      data = leaf_chi_subset)
summary(leafNarea_lmer) # N = 1,561
Anova(leafNarea_lmer)

### assign treatment group labels
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '0' & leaf_chi_subset$Ktrt_fac == '0'] <- '-P, -K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '1' & leaf_chi_subset$Ktrt_fac == '0'] <- '+P, -K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '0' & leaf_chi_subset$Ktrt_fac == '1'] <- '-P, +K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '1' & leaf_chi_subset$Ktrt_fac == '1'] <- '+P, +K'
    
leafNarea_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                                NPgroup = c('-P, -K', '-P, -K', '+P, -K', '+P, -K', 
                                            '-P, +K', '-P, +K', '+P, +K', '+P, +K'),
                                Ntrt_fac = c(0, 1, 0, 1, 0, 1, 0, 1),
                                y = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5),
                                group = c(cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[1, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[2, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[3, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[4, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[5, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[6, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[7, 9],
                                          cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[8, 9]))
leafNarea_letters$Ntrt_fac <- as.factor(leafNarea_letters$Ntrt_fac)
leafNarea_letters$letter[leafNarea_letters$group == " 1  "] <- "a"
leafNarea_letters$letter[leafNarea_letters$group == " 12 "] <- "ab"
leafNarea_letters$letter[leafNarea_letters$group == "  2 "] <- "b"
leafNarea_letters$letter[leafNarea_letters$group == "   3"] <- "c"

(narea_plot <- ggplot(data = leaf_chi_subset, 
                         aes(x = PKgroup, y = log(narea), fill = Ntrt_fac)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 20, colour = 'black'),
        axis.text.y = element_text(size = 20, colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_text(data = leafNarea_letters, aes(x = x, y = y, label = letter), size = 6) +
    scale_fill_manual(values = c("gray40", "burlywood1"), labels = c("Ambient", "Added N")) +
    labs(fill = "Soil N") +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab('P x K treatment'))

## find relative importance for each factor from model
relimp_leafn <- NULL
relimp_leafn$Factor <- c('Soil N', 'Soil P', 'Soil K+µ', 'Temperature', 'PAR', 'LMA', 'χ',
                         'N fixer', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafn$Importance <- as.numeric(as.character(c(calc.relip.mm(leafNarea_lmer)$lmg[1:9], 
                                                        sum(calc.relip.mm(leafNarea_lmer)$lmg[10:13]),
                                                        1 - sum(calc.relip.mm(leafNarea_lmer)$lmg))))

relimp_leafn_df <- as.data.frame(relimp_leafn)

tm <- treemapify(data = relimp_leafn_df,
                 area = "Importance", start = "topleft")
tm$x <- (tm$xmax + tm$xmin) / 2
tm$y <- (tm$ymax + tm$ymin) / 2

narea_tm <- full_join(relimp_leafn_df, tm, by = "Factor")
narea_tm$name <- c('Soil~N', 'Soil~P', 'Soil~K[+µ]', 'italic(T)[g]', 'italic(I)[g]',
                   'italic(M)[area]', 'χ', 'N~fixer', 
                   'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')

(narea_treemap <- ggplot(narea_tm, 
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                             label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "right",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    scale_fill_gradient(low = "lightcyan", high = "lightcyan4") +
    geom_text(data = filter(narea_tm, Factor == 'LMA' | Factor == 'PAR'), 
              aes(x = x, y = y), parse = T, size = 14) +
    geom_text(data = filter(narea_tm, Factor == 'χ'), 
              aes(x = x, y = y), parse = T, size = 10, family = 'Times') +
    geom_text(data = filter(narea_tm, Factor == 'Unexplained' | Factor == 'Temperature' | 
                              Factor == 'N fixer' | Factor == 'C3/C4'), 
              aes(x = x, y = y), parse = T, size = 7) +
    geom_text(data = filter(narea_tm, Factor == 'Soil N'), 
              aes(x = x, y = y), parse = T, size = 4) +
    ggrepel::geom_text_repel(data = filter(narea_tm, Factor == 'Soil Interactions' | Factor == 'Soil P' |
                                             Factor == 'Soil K+µ'), 
                             aes(x = x, y = y), parse = T, size = 4, 
                             direction = "y", xlim = c(1.01, NA)) +
    scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)))

#### impacts of nitrogen demand and nitrogen availability on Narea ####
### calculate vcmax, jmax, and vpmax with known chi for C3 plants
leaf_chi_subset_c3 <- subset(leaf_chi_subset, photosynthetic_pathway == 'C3')
gas_exchange_pred_c3 <- calc_optimal_vcmax(pathway = "C3",
                                           tg_c = leaf_chi_subset_c3$tmp, 
                                           paro = leaf_chi_subset_c3$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c3$vpd, 
                                           z = leaf_chi_subset_c3$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c3$chi,
                                           lma = leaf_chi_subset_c3$lma)

### calculate vcmax, jmax, and vpmax with known chi for C4 plants
leaf_chi_subset_c4 <- subset(leaf_chi_subset, photosynthetic_pathway == 'C4')
gas_exchange_pred_c4 <- calc_optimal_vcmax(pathway = "C4",
                                           tg_c = leaf_chi_subset_c4$tmp, 
                                           paro = leaf_chi_subset_c4$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c4$vpd, 
                                           z = leaf_chi_subset_c4$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c4$chi,
                                           lma = leaf_chi_subset_c4$lma)

## add C3 model results to leaf dataset
npred_c3 <- bind_cols(leaf_chi_subset_c3, gas_exchange_pred_c3[ ,39:51])
npred_c3$model_lma <- gas_exchange_pred_c3$lma

## add C4 model results to leaf dataset
npred_c4 <- bind_cols(leaf_chi_subset_c4, gas_exchange_pred_c4[ ,39:51])
npred_c4$model_lma <- gas_exchange_pred_c4$lma

# combine C3 and C4 subsets
leaf_chi_subset_npred <- bind_rows(npred_c3, npred_c4)

leaf_chi_subset_npred$lognphoto <- log(leaf_chi_subset_npred$nphoto)
leaf_chi_subset_npred$lognstructure <- log(leaf_chi_subset_npred$nstructure)

### fit linear mixed effects model
npred_soil_lmer <- lmer(log(narea) ~ lognphoto + lognstructure +
                          Ntrt_fac * Ptrt_fac * Ktrt_fac +
                          Nfix + photosynthetic_pathway + 
                          (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),
                        data = leaf_chi_subset_npred)
summary(npred_soil_lmer) # N = 1,548
Anova(npred_soil_lmer)

### make figures
## find slope and intercept from mixed effects model
nphoto_slope <- summary(emtrends(npred_soil_lmer, ~lognphoto, var = "lognphoto"))[1, 2]
nphoto_intercept_lowN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[1, 2]
nphoto_intercept_highN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[2, 2]
lognphoto_seq <- seq(min(leaf_chi_subset_npred$lognphoto, na.rm = T), 
                     max(leaf_chi_subset_npred$lognphoto, na.rm = T), 0.01)
nphoto_trend_lowN <- nphoto_intercept_lowN + lognphoto_seq * nphoto_slope
nphoto_trend_highN <- nphoto_intercept_highN + lognphoto_seq * nphoto_slope
nphoto_trend <- as.data.frame(cbind(lognphoto_seq, nphoto_trend_lowN, nphoto_trend_highN))

(npred_photo_plot <- ggplot(data = leaf_chi_subset_npred, 
                            aes(x = lognphoto, y = log(narea), color = Ntrt_fac)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    geom_point(shape = 16, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
    labs(color = 'Soil N') +
    geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_lowN), 
              col = 'black', lwd = 2, alpha = 0.8) +
    geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_highN), 
              col = 'burlywood1', lwd = 3, alpha = 0.8) +
    scale_y_continuous(limits = c(-2.5, 5.5)) +
    scale_x_continuous(limits = c(-2, 0)) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('N')['photo'])))
  
## find slope and intercept from mixed effects model
nstructure_slope <- summary(emtrends(npred_soil_lmer, ~lognstructure, var = "lognstructure"))[1, 2]
nstructure_intercept_lowN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[1, 2]
nstructure_intercept_highN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[2, 2]
lognstructure_seq <- seq(min(leaf_chi_subset_npred$lognstructure, na.rm = T), 
                         max(leaf_chi_subset_npred$lognstructure, na.rm = T), 0.1)
nstructure_trend_lowN <- nstructure_intercept_lowN + lognstructure_seq * nstructure_slope
nstructure_trend_highN <- nstructure_intercept_highN + lognstructure_seq * nstructure_slope
nstructure_trend <- as.data.frame(cbind(lognstructure_seq, nstructure_trend_lowN, nstructure_trend_highN))

(npred_structure_plot <- ggplot(data = leaf_chi_subset_npred, 
                              aes(x = lognstructure, y = log(narea), color = Ntrt_fac)) +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key=element_blank(),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    geom_point(shape = 16, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
    labs(color = 'Soil N') +
    geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_lowN), 
              col = 'black', lwd = 3, alpha = 0.8) +
    geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_highN), 
              col = 'burlywood1', lwd = 3, alpha = 0.8) +
    scale_y_continuous(limits = c(-2.5, 5.5)) +
    scale_x_continuous(limits = c(-5, 2.75)) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('N')['structure'])))

(npred_plot <- npred_structure_plot + npred_photo_plot +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 18)))

## find relative importance for each factor from model
relimp_leafn_pred <- NULL
relimp_leafn_pred$Factor <- c('N photo', 'N structure', 'Soil N', 'Soil P', 'Soil K+µ', 
                              'N fixer', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafn_pred$Importance <- as.numeric(as.character(c(calc.relip.mm(npred_soil_lmer)$lmg[1:7], 
                                                     sum(calc.relip.mm(npred_soil_lmer)$lmg[8:11]),
                                                     1 - sum(calc.relip.mm(npred_soil_lmer)$lmg))))
relimp_leafn_pred_df <- as.data.frame(relimp_leafn_pred)

sum(calc.relip.mm(npred_soil_lmer)$lmg[c(3:5, 8:11)])

tm_pred <- treemapify(data = relimp_leafn_pred_df,
                 area = "Importance", start = "topleft")
tm_pred$x <- (tm_pred$xmax + tm_pred$xmin) / 2
tm_pred$y <- (tm_pred$ymax + tm_pred$ymin) / 2

narea_tm_pred <- full_join(relimp_leafn_pred_df, tm_pred, by = "Factor")
narea_tm_pred$name <- c('italic(N)[photo]', 'italic(N)[structure]', 'Soil~N', 'Soil~P', 'Soil~K[+µ]', 
                        'N~fixer', 'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')

(narea_pred_treemap <- ggplot(narea_tm_pred, 
                              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                                  label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "right",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    scale_fill_gradient(low = "lightcyan", high = "lightcyan4") +
    geom_text(data = filter(narea_tm_pred, Factor == 'Unexplained'),
              aes(x = x, y = y), parse = T, size = 14) +
    geom_text(data = filter(narea_tm_pred, Factor == 'N structure' | Factor == 'N photo' | 
                              Factor == 'Soil N'),
              aes(x = x, y = y), parse = T, size = 10) +
    geom_text(data = filter(narea_tm_pred, Factor == 'Soil K+µ' | Factor == 'Soil P' | 
                              Factor == 'N fixer'),
              aes(x = x, y = y), parse = T, size = 8) +
    geom_text(data = filter(narea_tm_pred, Factor == 'Soil Interactions'),
              aes(x = x, y = y), parse = T, size = 7) +
    geom_text(data = filter(narea_tm_pred, Factor == 'C3/C4'),
              aes(x = x, y = y), parse = T, size = 6))

#### response of aboveground biomass to the soil treatments ####
### load data
leaf_site <- read.csv('../Data/site_traits.csv')

## turn treatment numbers into factors
leaf_site$Ntrt_fac <- as.factor(leaf_site$Ntrt)
leaf_site$Ptrt_fac <- as.factor(leaf_site$Ptrt)
leaf_site$Ktrt_fac <- as.factor(leaf_site$Ktrt)
leaf_site$block_fac <- as.factor(leaf_site$block)

leaf_site_sub = subset(leaf_site, site_code == "bldr.us" | site_code == "bogong.au" | 
                         site_code == "burrawan.au" | site_code == "cbgb.us" | site_code == "comp.pt" | 
                         site_code == "cowi.ca" | site_code == "elliot.us" | 
                         site_code == "gilb.za" | site_code == "hopl.us" | site_code == "lancaster.uk" | 
                         site_code == "mcla.us" | site_code == "mtca.au" | 
                         site_code == "sage.us" | site_code == "sgs.us" | 
                         site_code == "sier.us" | site_code == "smith.us" | 
                         site_code == "summ.za" | site_code == "unc.us" | site_code == "valm.ch")

### linear mixed effects model for mean live mass (AGB)
live_mass_lmer <- lmer(log(live_mass_mean) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac +
                        (1|site_code) + (1|site_code:block), 
                      data = leaf_site_sub)
summary(live_mass_lmer) # N = 763
Anova(live_mass_lmer)

### assign treatment group labels
leaf_site$PKgroup[leaf_site$Ptrt_fac == '0' & leaf_site$Ktrt_fac == '0'] <- '-P, -K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '1' & leaf_site$Ktrt_fac == '0'] <- '+P, -K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '0' & leaf_site$Ktrt_fac == '1'] <- '-P, +K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '1' & leaf_site$Ktrt_fac == '1'] <- '+P, +K'

### make figures
live_mass_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                                PKgroup = c('-P, +K', '-P, -K', '-P, -K', '+P, +K',
                                            '+P, -K', '-P, +K', '+P, -K', '+P, +K'),
                                Ntrt_fac = c(0, 0, 1, 0, 0, 1, 1, 1),
                                y = c(7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8), 
                                group = c(cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[1, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[2, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[3, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[4, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[5, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[6, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[7, 9],
                                          cld(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))[8, 9]))
live_mass_letters$Ntrt_fac <- as.factor(live_mass_letters$Ntrt_fac)
live_mass_letters$letter[live_mass_letters$group == " 1  "] <- "a"
live_mass_letters$letter[live_mass_letters$group == " 12 "] <- "ab"
live_mass_letters$letter[live_mass_letters$group == "  2 "] <- "b"
live_mass_letters$letter[live_mass_letters$group == "  23"] <- "bc"
live_mass_letters$letter[live_mass_letters$group == "   3"] <- "c"

(live_mass_plot <- ggplot(data = leaf_site, 
                          aes(x = PKgroup, y = log(live_mass_mean), fill = Ntrt_fac)) +
    theme(legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    scale_fill_manual(values = c("gray40", "burlywood1"), labels = c("Ambient", "Added N")) +
    scale_y_continuous(limits = c(3.5, 8)) +
    scale_x_discrete(limits = c('-P, -K', '-P, +K', '+P, -K', '+P, +K')) +
    labs(fill = "Soil N") +
    xlab('P x K treatment') +
    ylab(expression('ln AGB (g)')))

#### effect of soil nitrogen addition in relation to community nitrogen demand on Narea ####
### calculate treatment type averages
leaf_site_N <- leaf_chi_subset_npred %>%
  group_by(site_code, Ntrt_fac, Ptrt_fac, Ktrt_fac, block_fac, fence, trt, Taxon) %>%
  summarise_at(vars(narea, spp_lai, chi, spp_live_mass, spp_mass_N, max_cover, lma, p_pet, vpd, tmp,
                    Ambient_PAR, Ground_PAR, par_rat),
               mean, na.rm = TRUE)

### subset by low N treatment
leaf_site_lowN <- subset(leaf_site_N, Ntrt_fac == '0')
nrow(leaf_site_lowN)

### subset by high N treatment
leaf_site_highN <- subset(leaf_site_N, Ntrt_fac == '1')
nrow(leaf_site_highN)

### combine low N and high N subsets
leaf_site_trt <- left_join(leaf_site_lowN, leaf_site_highN, 
                                by = c('site_code', 
                                       'block_fac', 
                                       'Ptrt_fac', 'Ktrt_fac',
                                       'fence',
                                       'Taxon'))

### calculate percent change from the high N plots to the low N plots
leaf_site_trt$delta_narea <- ((leaf_site_trt$narea.y - 
                                leaf_site_trt$narea.x) / leaf_site_trt$narea.x) * 100
leaf_site_trt$delta_lai <- ((leaf_site_trt$spp_lai.y - 
                              leaf_site_trt$spp_lai.x) / leaf_site_trt$spp_lai.x) * 100
leaf_site_trt$delta_live_mass <- ((leaf_site_trt$spp_live_mass.y - 
                                    leaf_site_trt$spp_live_mass.x) / leaf_site_trt$spp_live_mass.x) * 100
leaf_site_trt$delta_N_mass <- ((leaf_site_trt$spp_mass_N.y - 
                                    leaf_site_trt$spp_mass_N.x) / leaf_site_trt$spp_mass_N.x) * 100
leaf_site_trt$delta_chi <- ((leaf_site_trt$chi.y - 
                              leaf_site_trt$chi.x) / leaf_site_trt$chi.x) * 100
leaf_site_trt$delta_lma <- ((leaf_site_trt$lma.y - 
                              leaf_site_trt$lma.x) / leaf_site_trt$lma.x) * 100

nrow(subset(leaf_site_trt, delta_narea > -1000))

### find mean absolute deviation (Leys et al., 2013)
delta_narea_mad <- mad(leaf_site_trt$delta_narea, na.rm = T)
delta_lai_mad <- mad(leaf_site_trt$delta_lai, na.rm = T)
delta_lma_mad <- mad(leaf_site_trt$delta_lma, na.rm = T)
delta_live_mass_mad <- mad(leaf_site_trt$delta_live_mass, na.rm = T)
delta_N_mass_mad <- mad(leaf_site_trt$delta_N_mass, na.rm = T)
delta_chi_mad <- mad(leaf_site_trt$delta_chi, na.rm = T)

delta_live_mass_data <- subset(leaf_site_trt,
                             delta_narea < 3 * delta_narea_mad &
                               delta_narea > 3 * -delta_narea_mad &
                               delta_lma < 3 * delta_lma_mad &
                               delta_lma > 3 * -delta_lma_mad &
                               delta_live_mass < 3 * delta_live_mass_mad &
                               delta_live_mass > 3 * -delta_live_mass_mad)

### linear mixed effects model for delta Narea by delta live mass
delta_live_mass_lm <- lmer(delta_narea ~ 
                            delta_live_mass + 
                            Ptrt_fac * Ktrt_fac +
                            delta_lma +
                            delta_live_mass : delta_lma +
                            (1|Taxon) + (1|Taxon:site_code) +
                            (1|Taxon:block_fac:site_code), 
                        data = delta_live_mass_data)
summary(delta_live_mass_lm) # N = 328
Anova(delta_live_mass_lm)

### make figures
## dataset
delta_live_mass_plot_data <- delta_live_mass_data

## trendline information
delta_live_mass_plot_intercept_lowlma <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass, 
                                                         at = list(delta_live_mass = 0, delta_lma = -25)))[1, 2]
delta_live_mass_plot_slope_lowlma <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
                                                      at = list(delta_lma = -25)))[1, 2]
delta_live_mass_plot_intercept_midlma <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass, 
                                                         at = list(delta_live_mass = 0, delta_lma = 0)))[1, 2]
delta_live_mass_plot_slope_midlma <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
                                                      at = list(delta_lma = 0)))[1, 2]
delta_live_mass_plot_intercept_highlma <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass, 
                                                          at = list(delta_live_mass = 0, delta_lma = 25)))[1, 2]
delta_live_mass_plot_slope_highlma <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
                                                       at = list(delta_lma = 25)))[1, 2]

delta_live_mass_plot_trend_lowlma <- delta_live_mass_plot_slope_lowlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_lowlma

delta_live_mass_plot_trend_midlma <- delta_live_mass_plot_slope_midlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_midlma

delta_live_mass_plot_trend_highlma <- delta_live_mass_plot_slope_highlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_highlma

delta_live_mass_plot_trend_df <- data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass), 
                                                max(delta_live_mass_plot_data$delta_live_mass), 1),
                                            delta_live_mass_plot_trend_lowlma, 
                                            delta_live_mass_plot_trend_midlma, 
                                            delta_live_mass_plot_trend_highlma)
colnames(delta_live_mass_plot_trend_df) <- c('delta_live_mass', 
                                             'delta_narea_lowlma', 
                                             'delta_narea_midlma', 
                                             'delta_narea_highlma')

(delta_live_mass_plot <- ggplot(data = delta_live_mass_plot_data, 
                                aes(x = delta_live_mass, y = delta_narea, fill = delta_lma, size = delta_lma)) +
    theme(axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    geom_point(shape = 21, colour = 'black', stroke = 0.5, alpha = 1) +
    scale_size_continuous(range = c(0.5, 3)) +
    scale_fill_gradient(low = 'grey80', high = 'grey0') +
    geom_line(data = delta_live_mass_plot_trend_df, 
              aes(x = delta_live_mass, y = delta_narea_lowlma, fill = NULL), 
              size = 4, colour = 'grey60', alpha = 1, lty = 2) +
    geom_line(data = delta_live_mass_plot_trend_df, 
              aes(x = delta_live_mass, y = delta_narea_midlma, fill = NULL), 
              size = 4, colour = 'grey40', alpha = 1, lty = 2) +
    geom_line(data = delta_live_mass_plot_trend_df, 
              aes(x = delta_live_mass, y = delta_narea_highlma, fill = NULL), 
              size = 4, colour = 'grey20', alpha = 1, lty = 1) +
    labs(fill = expression('∆' * italic('M')['area'] * ' (%)')) +
    guides(size = "none") +
    ylab(expression('∆' * italic('N')['area'] * ' (%)')) +
    xlab(expression('∆' *'AGB' * ' (%)')))

