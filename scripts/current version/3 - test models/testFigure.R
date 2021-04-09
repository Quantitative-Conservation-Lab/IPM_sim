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
