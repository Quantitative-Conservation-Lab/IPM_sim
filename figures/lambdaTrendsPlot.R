# AEB NOTE 
# this needs to be modified

# 9 panel plot
# rows - low, med, high lambda
# columns - low, med, high detection

#or 3? column is model type?

library(tidybayes)
library(ggplot2)
library(wesanderson)

pal <- rev(wes_palette("Zissou1", 3, type = "continuous"))

# AEB - idea for version two of this plot
## row low
g1 <- ggplot(toplot1, aes(x = Year, y = value)) +
 stat_pointinterval(aes(color = simscenarios, fill = simscenarios), alpha = 0.5, .width = c(0.5, 0.95)) +
 geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
 geom_hline(yintercept = 1.00, linetype = "solid", color = "black") +
 facet_grid(.~simscenarios) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "none", # turn off legend
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_discrete(name="Year",breaks = seq(1, 14, by = 3)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))

## row med
g2 <- ggplot(toplot2, aes(x = Year, y = value)) +
  stat_pointinterval(aes(color = simscenarios, fill = simscenarios), alpha = 0.5, .width = c(0.5, 0.95)) +
  geom_hline(yintercept = 1.00, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black") +
  facet_grid(.~simscenarios) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "none", # turn off legend
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_discrete(name="Year",breaks = seq(1, 14, by = 3)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))

## row high
g3 <- ggplot(toplot3, aes(x = Year, y = value)) +
  stat_pointinterval(aes(color = simscenarios, fill = simscenarios), alpha = 0.5, .width = c(0.5, 0.95)) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black") +
  facet_grid(.~simscenarios) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  scale_x_discrete(name="Year", labels = seq(1, 14, by = 3), breaks = seq(1, 14, by = 3)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))

## combine all
library(patchwork)
library(cowplot)

pdf(here("figures", "lambda_plot1.pdf"), width = 12, height = 8)
plot_grid(g1, g2, g3, nrow = 3)
dev.off()


# AEB - old stuff ########
p1 <- ggplot(transform(toplot1,
                       model = factor(model, levels = c("Full IPM", "No nest data", "No mark recapture data", "Abundance data only")))) +
  #geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_line(aes(x = Year, y = high, color = detection), linetype = "solid", alpha = 0.5) +
  geom_line(aes(x = Year, y = low, color = detection), linetype = "solid", alpha = 0.5) +
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
  scale_x_continuous(name="Year", limits = c(1, 13), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name="", limits = c(0.8, 1.3), breaks = c(0.95, 1, 1.05), position = "left") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925+0.1, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
p1

p2 <- ggplot(toplot2) +
  #geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_line(aes(x = Year, y = high, color = detection), linetype = "solid", alpha = 0.5) +
  geom_line(aes(x = Year, y = low, color = detection), linetype = "solid", alpha = 0.5) +
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
  scale_x_continuous(name="Year", limits = c(1, 13), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = c(0.95, 1, 1.05), position = "left") +
  #scale_x_continuous(name="Year", limits = c(1, 13), breaks = 1:14) +
  #scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.1)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925+0.05, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
p2

p3 <- ggplot(toplot3) +
  #geom_hline(yintercept = 1, linetype = "solid", color = "grey40") +
  geom_line(aes(x = Year, y = high, color = detection), linetype = "solid", alpha = 0.5) +
  geom_line(aes(x = Year, y = low, color = detection), linetype = "solid", alpha = 0.5) +
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
  scale_x_continuous(name="Year", limits = c(1, 13), breaks = seq(1, 14, by = 3)) +
  scale_y_continuous(name="", limits = c(0.8, 1.3), breaks = c(0.95, 1, 1.05), position = "left") +
  #scale_x_continuous(name="Year", limits = c(1, 13), breaks = 1:14) +
  #scale_y_continuous(name=expression(paste("Estimated ", lambda)), limits = c(0.8, 1.3), breaks = seq(0.8, 1.3, 0.1)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) +
  coord_cartesian(clip = "off") +
  #annotate(geom = "text", x=14.5, y=0.65, label="Population\ndeclining", size = 3) +
  #annotate(geom = "text", x=14.5, y=1.25, label="Population\nincreasing", size = 3) +
  annotate(geom = "text", x=14.5, y=0.925, label=expression(paste("True ",lambda)), size = 3) +
  scale_color_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Data\nquality", labels = c("Low", "Medium", "High"))
p3

library(patchwork)
library(cowplot)

pdf("IPMplotTrevor.pdf", width = 12, height = 8)
plot_grid(p1, p2, p3, nrow = 3)
dev.off()
