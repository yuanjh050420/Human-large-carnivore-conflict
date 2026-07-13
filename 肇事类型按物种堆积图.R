library(ggplot2)
library(dplyr)
library(tidyr)

# Create wildlife conflict dataframe
data <- data.frame(
  Species = c("Asiatic black bear", "Grey wolf", "Tiger", "Brown bear", 
              "Snow leopard", "Leopard", "Eurasian lynx",  
              "Asiatic golden cat","Clouded leopard", "Dhole"),
  Number_of_events = c(187, 79, 78, 64, 44, 44, 16, 4, 2, 2),
  Attack_livestock = c(23, 66, 58, 18, 29, 30, 6, 2, 0, 0),
  Enter_human_settlements = c(34, 12, 15, 17, 14, 13, 10, 2, 2, 2),
  Attack_human = c(75, 1, 5, 9, 0, 1, 0, 0, 0, 0),
  Cause_property_damage = c(55, 0, 0, 20, 1, 0, 0, 0, 0, 0)
)

# Convert to long format for ggplot
data_long <- data %>%
  select(Species, Attack_livestock, Enter_human_settlements, 
         Attack_human, Cause_property_damage) %>%
  pivot_longer(cols = -Species, 
               names_to = "Conflict_Type", 
               values_to = "Count") %>%
  mutate(Conflict_Type = factor(Conflict_Type, 
                                levels = c("Attack_livestock", "Attack_human", 
                                           "Enter_human_settlements", "Cause_property_damage"),
                                labels = c("Livestock/dog attack", 
                                           "Human attack", 
                                           "Human settlement entry", 
                                           "Property damage")))

# Set species order (descending by total events)
species_order <- data$Species[order(data$Number_of_events, decreasing = TRUE)]
data_long$Species <- factor(data_long$Species, levels = species_order)

# Create stacked bar chart
p <- ggplot(data_long, aes(x = Species, y = Count, fill = Conflict_Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Number of events", fill = "Conflict type") +
  scale_x_discrete() +
  scale_fill_manual(values = c("Livestock/dog attack" = "#D55E00",
                               "Human attack" = "#FFB6C1",
                               "Human settlement entry" = "#0072B2",
                               "Property damage" = "#90EE90"),
                    breaks = c("Livestock/dog attack", "Human attack", 
                               "Human settlement entry", "Property damage")) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, hjust = 1, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 10, 10, 20)
  )

# Display the plot
print(p)

ggsave("C:/Users/yuanj/Desktop/人兽冲突/整体作图/肇事类型按物种堆积/wildlife_conflict_stacked_bar.png", p,
        width = 11, height = 6, units = 'in',
        dpi = 300, bg = "white")

