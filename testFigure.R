# aeb
# march 2, 2020
# example post processing
library(tidyverse)
library(here)
library(nimble)

# load data
scenarios <- readRDS(here("data", "scenarios.RDS"))
low.lam.combos <- readRDS(here("data", "low.lam.combos.RDS"))
med.lam.combos <- readRDS(here("data", "med.lam.combos.RDS"))
high.lam.combos <- readRDS(here("data", "high.lam.combos.RDS"))

# functions
source(here("scripts", "current version",
            "1 - simulating data", "IPM_sim_2.0function.R"))
source(here("scripts", "current version",
            "2 - models", "IPMNimble_v2.0.R"))
source(here("scripts", "current version",
            "4 - run models", "run_scenarios_helperFns.R"))

low.comb <- low.lam.combos[sample(1:5000, 1), 1:3]
med.comb <- med.lam.combos[sample(1:5000, 1), 1:3]
high.comb <- high.lam.combos[sample(1:5000, 1), 1:3]

lowpopTraj <- simPopTrajectory(n.years=15,
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150),
                               phi.1=as.numeric(low.comb[2]),
                               phi.ad=as.numeric(low.comb[3]),
                               f=as.numeric(low.comb[1]))

medpopTraj <- simPopTrajectory(n.years=15,
                               n.data.types=c(0.25,0.25,0.25),
                               age.init=c(150,150),
                               phi.1=as.numeric(med.comb[2]),
                               phi.ad=as.numeric(med.comb[3]),
                               f=as.numeric(med.comb[1]))

highpopTraj <- simPopTrajectory(n.years=15,
                                n.data.types=c(0.25,0.25,0.25),
                                age.init=c(150,150),
                                phi.1=as.numeric(high.comb[2]),
                                phi.ad=as.numeric(high.comb[3]),
                                f=as.numeric(high.comb[1]))

# simulate data
detect.l <- 0.3
detect.m <- 0.5
detect.h <- 0.8

detect <- c(detect.l, detect.m, detect.h)

nb <- 100000#0 #burn-in
ni <- nb + nb #total iterations
nt <- 10  #thin
nc <- 3  #chains

for (d in 1:3) {
  print(paste("detection level ", detect[d]))
  lowpopDat <- simData (indfates = lowpopTraj$indfates,
                        n.years = 15,
                        n.data.types = c(0.25,0.25,0.25),
                        ADonly = T,
                        p.1 = detect[d],
                        p.ad = detect[d],
                        p.count = detect[d],
                        p.prod = detect[d],
                        BinMod = T,
                        n.sam = 3,
                        sig = 0,
                        productivity = T)
  medpopDat <- simData (indfates = medpopTraj$indfates,
                        n.years = 15,
                        n.data.types = c(0.25,0.25,0.25),
                        ADonly = T,
                        p.1 = detect[d],
                        p.ad = detect[d],
                        p.count = detect[d],
                        p.prod = detect[d],
                        BinMod = T,
                        n.sam = 3,
                        sig = 0,
                        productivity = T)
  highpopDat <- simData (indfates = highpopTraj$indfates,
                         n.years = 15,
                         n.data.types = c(0.25,0.25,0.25),
                         ADonly = T,
                         p.1 = detect[d],
                         p.ad = detect[d],
                         p.count = detect[d],
                         p.prod = detect[d],
                         BinMod = T,
                         n.sam = 3,
                         sig = 0,
                         productivity = T)
  for (m in 1:4) {
    print(paste("model ", m))
    if(m == 1) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runIPMmod(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 2) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runnonests(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 3) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runnomr(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    } else if (m == 4) {
      popDat <- highpopDat
      popTraj <- highpopTraj
      comb <- high.comb
      highout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- medpopDat
      popTraj <- medpopTraj
      comb <- med.comb
      medout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))

      popDat <- lowpopDat
      popTraj <- lowpopTraj
      comb <- low.comb
      lowout <- runabundonly(nb = nb, ni = ni, nt = nt, nc = nc, popDat, popTraj, comb, detect = rep(detect[d], 3))
    }
    assign(paste("highout-",d,"-",m, sep = ""), highout)
    assign(paste("medout-",d,"-",m, sep = ""), medout)
    assign(paste("lowout-",d,"-",m, sep = ""), lowout)

    saveRDS(highout, paste("highout-",d,"-",m, ".RDS", sep = ""))
    saveRDS(medout, paste("medout-",d,"-",m,  ".RDS", sep = ""))
    saveRDS(lowout, paste("lowout-",d,"-",m,  ".RDS", sep = ""))
  }
}

#### Processing #######

# 12 panel plot
  # rows - low, med, high out
  # columns - model type
  # on each plot - 3 lines representing detection type

# row 1
# declining population

row1 <- bind_rows(as.data.frame(as.matrix(`highout-1-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-1-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-1-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-1-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-2-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-2-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-2-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-2-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-3-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-3-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-3-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`highout-3-4`)[1:1000, ])
                  ) %>%
  select(contains("lambda"))

for (i in 1:14) {
  row1 <- row1 %>%
    mutate("geomean.{i}" :=  NA_real_)
}

for(i in 1:dim(row1)[1]) {
  for(j in 1:((ncol(row1) - 2 )/2)) {
    print(paste("row", i, "col", j))
    row1[i, (ncol(row1)/2 + j)] <- exp(mean(unlist(log(row1[i, 1:j]))))
  }
}

row1 <- cbind(row1,
              model = rep(rep(1:4, each = 1000), times = 3),
              detection = rep(1:3, each = 1000*4)) %>%
  as.data.frame()

toplot1 <- row1 %>%
  select(contains("geomean"), model, detection) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(model, detection, Year) %>%
  summarise(low = quantile(value, 0.025),
            med = quantile(value, 0.5),
            high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(model = as.factor(model),
         detection = as.factor(detection)) %>%
  mutate(model = case_when(
    model == 1 ~ "Full IPM",
    model == 2 ~ "No nest data",
    model == 3 ~ "No mark recapture data",
    model == 4 ~ "Abundance data only"
  ))
levels(toplot1$model) = c("Full IPM", "No nest data", "No mark recapture data", "Abundance data only")

library(ggplot2)
library(gtable)
library(RColorBrewer)
library(wesanderson)

pal <- rev(wes_palette("Zissou1", 3, type = "continuous"))

#pdf("IPMplotTrevor.pdf", width = 12, height = 8)
p1 <- ggplot(transform(toplot1,
                       model = factor(model, levels = c("Full IPM", "No nest data", "No mark recapture data", "Abundance data only")))) +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_ribbon(aes(x = Year, ymin = low, ymax = high, fill = detection), alpha = 0.5) +
  geom_line(aes(x = Year, y = med, color = detection), size = 1.05) +
  facet_grid(.~model) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_continuous(name="Year", limits = c(1, 16), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name="", limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.2)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925+0.1, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
p1
#dev.off()

# row 2
# steady population

row2 <- bind_rows(as.data.frame(as.matrix(`medout-1-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-1-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-1-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-1-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-2-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-2-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-2-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-2-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-3-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-3-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-3-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`medout-3-4`)[1:1000, ])
) %>%
  select(contains("lambda"))

for (i in 1:14) {
  row2 <- row2 %>%
    mutate("geomean.{i}" :=  NA_real_)
}

for(i in 1:dim(row2)[1]) {
  for(j in 1:((ncol(row2) - 2 )/2)) {
    print(paste("row", i, "col", j))
    row2[i, (ncol(row2)/2 + j)] <- exp(mean(unlist(log(row2[i, 1:j]))))
  }
}

row2 <- cbind(row2,
              model = rep(rep(1:4, each = 1000), times = 3),
              detection = rep(1:3, each = 1000*4)) %>%
  as.data.frame()

toplot2 <- row2 %>%
  select(contains("geomean"), model, detection) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(model, detection, Year) %>%
  summarise(low = quantile(value, 0.025),
            med = quantile(value, 0.5),
            high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(model = as.factor(model),
         detection = as.factor(detection)) %>%
  mutate(high = if_else(high > 1.3, 1.299, high))


library(ggplot2)
library(gtable)
library(RColorBrewer)
#pdf("IPMplotTrevor.pdf", width = 12, height = 8)
p2 <- ggplot(toplot2) +
  #geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_ribbon(aes(x = Year, ymin = low, ymax = high, fill = detection), alpha = 0.5) +
  geom_line(aes(x = Year, y = med, color = detection), size = 1.05) +
  facet_grid(.~model) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  #geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_continuous(name="Year", limits = c(1, 16), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.2)) +
  #scale_x_continuous(name="Year", limits = c(1, 16), breaks = 1:14) +
  #scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.1)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925+0.05, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
#p2

# row 3
# increasing population
row3 <- bind_rows(as.data.frame(as.matrix(`lowout-1-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-1-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-1-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-1-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-2-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-2-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-2-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-2-4`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-3-1`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-3-2`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-3-3`)[1:1000, ]),
                  as.data.frame(as.matrix(`lowout-3-4`)[1:1000, ])
) %>%
  select(contains("lambda"))

for (i in 1:14) {
  row3 <- row3 %>%
    mutate("geomean.{i}" :=  NA_real_)
}

for(i in 1:dim(row3)[1]) {
  for(j in 1:((ncol(row3) - 2 )/2)) {
    print(paste("row", i, "col", j))
    row3[i, (ncol(row3)/2 + j)] <- exp(mean(unlist(log(row3[i, 1:j]))))
  }
}

row3 <- cbind(row3,
              model = rep(rep(1:4, each = 1000), times = 3),
              detection = rep(1:3, each = 1000*4)) %>%
  as.data.frame()

toplot3 <- row3 %>%
  select(contains("geomean"), model, detection) %>%
  #group_by(model, detection)
  pivot_longer(cols = starts_with("geomean"), names_to = "Year") %>%
  filter(!is.na(value)) %>%
  mutate(Year = str_remove(Year, "geomean\\.")) %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(model, detection, Year) %>%
  summarise(low = quantile(value, 0.025),
            med = quantile(value, 0.5),
            high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(model = as.factor(model),
         detection = as.factor(detection))

library(ggplot2)
library(gtable)
library(RColorBrewer)

p3 <- ggplot(toplot3) +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_ribbon(aes(x = Year, ymin = low, ymax = high, fill = detection), alpha = 0.5) +
  geom_line(aes(x = Year, y = med, color = detection), size = 1.05) +
  facet_grid(.~model) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  #geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        legend.position = "none",
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_continuous(name="Year", limits = c(1, 16), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name="", limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.2)) +
  #scale_x_continuous(name="Year", limits = c(1, 16), breaks = 1:14) +
  #scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.1)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
#p3

library(patchwork)
library(cowplot)

pdf("IPMplotTrevor.pdf", width = 12, height = 8)
p1 / p2 / p3
dev.off()

# TODO
# fix figures
  # labels
  # legends
  # axes
  # text sizes
# figure out why the estimates are wrong :(


# Note for future ############
# in the future we should thin these by a lot (like down to 3000)8

# test1 <- as.data.frame(as.matrix(test1))
# test2 <- as.data.frame(as.matrix(test2))
# test3 <- as.data.frame(as.matrix(test3))
# test4 <- as.data.frame(as.matrix(test4))
# test5 <- as.data.frame(as.matrix(test5))
# test6 <- as.data.frame(as.matrix(test6))
# test7 <- as.data.frame(as.matrix(test7))
# test8 <- as.data.frame(as.matrix(test8))
# test9 <- as.data.frame(as.matrix(test9))
#
# out1 <-test1 %>%
#   select(contains("lambda"))
# out2 <- test2 %>%
#   select(contains("lambda"))
# out3 <- test3 %>%
#   select(contains("lambda"))
# out4 <- test4 %>%
#   select(contains("lambda"))
# out5 <- test5 %>%
#   select(contains("lambda"))
# out6 <- test6 %>%
#   select(contains("lambda"))
# out7 <- test7 %>%
#   select(contains("lambda"))
# out8 <- test8 %>%
#   select(contains("lambda"))
# out9 <- test9 %>%
#   select(contains("lambda"))
#
# #apply(tmp, 1, geoMean, na.rm = TRUE)
#
# # can do this easier in the tidyverse
#
# for (i in 1:ncol(out1)) {
#   out1 <- out1 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out2 <- out2 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out3 <- out3 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out4 <- out4 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out5 <- out5 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out6 <- out6 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out7 <- out7 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out8 <- out8 %>%
#     mutate("geomean.{i}" :=  NA_real_)
#   out9 <- out9 %>%
#     mutate("geomean.{i}" :=  NA_real_)
# }
#
# for(i in 1:1000) {
#   for(j in 1:(ncol(out1)/2)) {
#     print(paste("row", i, "col", j))
#     out1[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out1[i, 1:j]))))
#
#     out2[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out2[i, 1:j]))))
#
#     out3[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out3[i, 1:j]))))
#
#     out4[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out4[i, 1:j]))))
#
#     out5[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out5[i, 1:j]))))
#
#     out6[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out6[i, 1:j]))))
#
#     out7[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out7[i, 1:j]))))
#
#     out8[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out8[i, 1:j]))))
#
#     out9[i, (ncol(out1)/2 + j)] <- exp(mean(unlist(log(out9[i, 1:j]))))
#   }
# }
# beep(sound = 2)
#
# # sweet, now what percentage of the time are we estimating lambda < 1 in each year
#
# # truth is declining
# # proportion of time the geometric mean is less than 1
#
# # 1, 2, 3 should get grouped
# # 4, 5, 6 should get grouped
# # 7, 8, 9 should get grouped
#
#
# # what are the medians and quantiles by group
#
# toplot1 <- bind_rows(out1, out2, out3) %>%
#   select(contains("geomean")) %>%
#   pivot_longer(everything()) %>%
#   filter(!is.na(value)) %>%
#   mutate(name = str_remove(name, "geomean\\.")) %>%
#   mutate(name = as.numeric(name)) %>%
#   group_by(name) %>%
#   summarise(low = quantile(value, 0.025),
#             med = quantile(value, 0.5),
#             high = quantile(value, 0.975)) %>%
#   ungroup() %>%
#   bind_cols(group = rep(1, 14))
#
# toplot2 <- bind_rows(out4, out5, out6) %>%
#   select(contains("geomean")) %>%
#   pivot_longer(everything()) %>%
#   filter(!is.na(value)) %>%
#   mutate(name = str_remove(name, "geomean\\.")) %>%
#   mutate(name = as.numeric(name)) %>%
#   group_by(name) %>%
#   summarise(low = quantile(value, 0.025),
#             med = quantile(value, 0.5),
#             high = quantile(value, 0.975)) %>%
#   ungroup() %>%
#   bind_cols(group = rep(2, 14))
#
# toplot3 <- bind_rows(out7, out8, out9) %>%
#   select(contains("geomean")) %>%
#   pivot_longer(everything()) %>%
#   filter(!is.na(value)) %>%
#   mutate(name = str_remove(name, "geomean\\.")) %>%
#   mutate(name = as.numeric(name)) %>%
#   group_by(name) %>%
#   summarise(low = quantile(value, 0.025),
#             med = quantile(value, 0.5),
#             high = quantile(value, 0.975)) %>%
#   ungroup() %>%
#   bind_cols(group = rep(3, 14))
#
# dat <- bind_rows(toplot1, toplot2, toplot3) %>%
#   mutate(group = as.factor(group)) %>%
#   mutate(name = as.numeric(name))
#
# # PLOTTING ########
#
# library(ggplot2)
# library(gtable)
# library(RColorBrewer)
# pdf("IPMplotTrevor.pdf", width = 12, height = 8)
# p1 <- ggplot(dat) +
#   geom_ribbon(aes(x = name, ymin = low, ymax = high, fill = group), alpha = 0.5) +
#   geom_line(aes(x = name, y = med, color = group), size = 1.05) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
#   geom_hline(yintercept = 1, linetype = "solid", color = "black") +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),
#         legend.position = "top",
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         plot.title.position = "plot",
#         axis.title=element_text(size=14)) +
#   scale_x_continuous(name="Year", limits = c(1, 16), breaks = 1:14) +
#   scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.1)) +
#   theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")) +
#   coord_cartesian(clip = "off") +
#   annotate(geom = "text", x=15.5, y=0.65, label="Population\ndeclining", size = 5) +
#   annotate(geom = "text", x=15.5, y=1.25, label="Population\nincreasing", size = 5) +
#   annotate(geom = "text", x=15.5, y=0.925, label=expression(paste("True ",lambda)), size = 5) +
#   scale_color_brewer(palette = "RdYlBu", name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
#   scale_fill_brewer(palette = "RdYlBu", name = "Data\nquality", labels = c("Low", "Medium", "High"))
# p1
# dev.off()
#
