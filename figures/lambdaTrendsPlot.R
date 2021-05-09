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

## simplified version of above plot
df1 <- toplot1 %>%
  mutate(Lambda = c("Decreasing"),
         simscenarios2 = ifelse(simscenarios == 1, "low", # check this
                                ifelse(simscenarios == 2, "med",
                                       ifelse(simscenarios == 3, "high", "")))) %>%
  slice_sample(prop = .1) # remove due to memory issues
df2 <- toplot2 %>%
  mutate(Lambda = c("Stable"),
         simscenarios2 = ifelse(simscenarios == 4, "low", # check this
                                ifelse(simscenarios == 5, "med",
                                       ifelse(simscenarios == 6, "high", "")))) %>%
  slice_sample(prop = .1) # remove due to memory issues
df3 <- toplot3 %>%
  mutate(Lambda = c("Increasing"),
         simscenarios2 = ifelse(simscenarios == 7, "low", # check this
                                ifelse(simscenarios == 8, "med",
                                       ifelse(simscenarios == 9, "high", "")))) %>%
  slice_sample(prop = .1) # remove due to memory issues
df <- rbind(df1, df2, df3)
# # change factor levels
# df$simscenarios <- factor(df$simscenarios, levels = c(1, 2, 3),
#                           ordered = TRUE, labels = c(expression(paste("low ", italic(p))), 
#                                                      expression(paste("med ", italic(p))), 
#                                                      expression(paste("high ", italic(p)))))
df$Lambda <- factor(df$Lambda, levels = c("Decreasing", "Stable", "Increasing"),
                                        ordered = TRUE, labels = c(expression(paste("Decreasing ", lambda)), 
                                                                   expression(paste("Stable ", lambda)), 
                                                                   expression(paste("Increasing ", lambda))))


# plotting
# h1 <- ggplot(df, aes(x = Year, y = value)) +  
#   stat_pointinterval(aes(color = simscenarios2, fill = simscenarios2), 
#                      alpha = 0.5, .width = c(0.5, 0.95)) +
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
#   geom_hline(yintercept = 1.00, linetype = "solid", color = "black") +
#   geom_hline(yintercept = 1.05, linetype = "dashed", color = "black") +
#   facet_grid(.~Lambda) +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),
#         # axis.text.x = element_blank(),
#         # axis.title.x = element_blank(),
#         # strip.text = element_blank(),
#         legend.position = "bottom", 
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 12),
#         plot.title.position = "plot",
#         axis.title=element_text(size=12)) +
#   scale_x_discrete(name="Year", breaks = seq(1, 14, by = 3)) +
#   coord_cartesian(clip = "off") +
#   scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
#   scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High")) 
# h1
h2 <- ggplot(df) +  
 geom_violin(aes(x = Year, y = value, color = simscenarios2, fill = simscenarios2, group = Year)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black") +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black") +
  facet_grid(.~Lambda) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        # strip.text = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.title=element_text(size=12)) +
  # scale_x_discrete(name="Year", breaks = seq(1, 14, by = 3)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High")) 
h2

# create plots for year 1, year 5, year 10 (need to re-factor lambda levels above)
dfyr1 <- df %>% filter(Year == 1)
dfyr5 <- df %>% filter(Year == 5)
dfyr10 <- df %>% filter(Year == 10)

## plotting
j1 <- ggplot(dfyr1) +  
  geom_violin(aes(x = simscenarios2, y = value, color = simscenarios2, fill = simscenarios2)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black", size = 0.25) +
  facet_wrap(.~Lambda, labeller = label_parsed, nrow = 3) +
  ylim(c(0.75, 1.25)) + ylab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        # strip.text = element_blank(),
        legend.position = "none", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        # plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5),
        axis.title=element_text(size=12)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High")) +
  ggtitle("Year 1")
j1
j2 <- ggplot(dfyr5) +  
  geom_violin(aes(x = simscenarios2, y = value, color = simscenarios2, fill = simscenarios2)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black", size = 0.25) +
  # facet_grid(.~Lambda, labeller = label_parsed) +
  facet_wrap(.~Lambda, labeller = label_parsed, nrow = 3) +
  ylim(c(0.75, 1.25)) + ylab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), # turn off y-axis
        axis.title.x = element_blank(),
        # strip.text = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        # plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5),
        axis.title=element_text(size=12)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High")) +
  ggtitle("Year 5")
j2
j3 <- ggplot(dfyr10) +  
  geom_violin(aes(x = simscenarios2, y = value, color = simscenarios2, fill = simscenarios2)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.00, linetype = "solid", color = "black", size = 0.25) +
  geom_hline(yintercept = 1.05, linetype = "dashed", color = "black", size = 0.25) +
  # facet_grid(.~Lambda, labeller = label_parsed) +
  facet_wrap(.~Lambda, labeller = label_parsed, nrow = 3) +
  ylim(c(0.75, 1.25)) + ylab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), # turn off y axis
        axis.title.x = element_blank(),
        # strip.text = element_blank(),
        legend.position = "none", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        # plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5),
        axis.title=element_text(size=12)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High"))  +
  scale_fill_manual(values = pal, name = "Detection\nlevel", labels = c("Low", "Medium", "High")) +
  ggtitle("Year 10")
j3

## combine
library(patchwork)
library(here)
all <- j1 + j2 + j3
all

## export out
ggsave(filename = here("figures", "lambda_plot4.pdf"), plot = all)

# ## try with cowplot
# library(cowplot)
# all2 <- align_plots(j1, j2, j3, align = "hv", axis = "tblr")
# all2
# 
# pdf("lambda_plot4.pdf", width = 8, height = 6)
# plot_grid(j1, j2, j3, ncol = 3)
# dev.off()

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
