# Load required packages
library(ggplot2)
library(dplyr)

dataset <- read.csv("clean.csv")

dataset_df <- dataset %>%
  group_by(TIME, CLEANING,SURFACE) %>%
  summarise(Mean_count = mean(COUNT, na.rm = TRUE),
            SE_count = sd(COUNT, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

a<-ggplot(dataset_df %>% filter(TIME>-.5), aes(x = TIME, y = Mean_count, color = CLEANING)) +
  geom_point(size=1.9) +
  geom_errorbar(aes(ymin = Mean_count - SE_count, ymax = Mean_count + SE_count), width = 0) +
  labs(x = "Time [hrs]", y = "Mean Number of Colonies",caption = "Mean data points \nerror bars represent standard error") +
  # facet_wrap(SURFACE~CLEANING,ncol = 1,nrow=3)+
  scale_color_brewer(palette = "Set2")+
  scale_y_continuous(trans = "log10",breaks = scales::breaks_pretty())+#+
  theme_minimal()
  # ggpubr::theme_classic2()


b<-ggplot(dataset_df %>% filter(TIME<0), aes(x = TIME, y = Mean_count, color = CLEANING)) +
  geom_point(size=1.9) +
  geom_errorbar(aes(ymin = Mean_count - SE_count, ymax = Mean_count + SE_count), width = 0) +
  labs(x = "Time [hrs]", y = "Mean Number of Colonies",caption = "Mean data points \nerror bars represent standard error") +
  # facet_wrap(SURFACE~CLEANING,ncol = 1,nrow=3)+
  scale_color_manual(values = c("red"))+
  scale_y_continuous(trans = "log10",breaks = scales::breaks_pretty())+#+
  theme_minimal()

ggpubr::ggarrange(b,a)


