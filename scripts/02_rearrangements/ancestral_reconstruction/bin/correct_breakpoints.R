library(dplyr)
library(ggplot2)

# Assume df is your data.frame with columns: chr, start, end, id, type, node, species

# Define node order (most recent first)

brp_file <- "combined_anc_breakpoints.bed"

brp.df <- read.csv(brp_file, sep = "\t", header = F)

colnames(brp.df) = c("chr", "start", "end", "brp_id", "type", "node", "species")

unique(brp.df$node)

# decides priority of node
node_levels <- c("N3", "N1") # adjust as needed

sink("tables/anc_breakpoint_counts_by_species_and_type_notes.txt", append = TRUE)
print("Summary of breakpoint counts before filtering by species and type:")
brp.df %>% group_by(species, type) %>% summarise(n = n())
sink()

View(brp.df)

df <- brp.df %>%
  mutate(node = factor(node, levels = node_levels)) %>%
  group_by(chr, start, end, type, species) %>%
  arrange(node, .by_group = TRUE) %>%
  summarise(
    event_parent_node = first(as.character(node)), # most recent node
    all_nodes = paste(sort(unique(as.character(node))), collapse = ","),
    all_brp_ids = paste(brp_id, collapse = ","),
    .groups = "drop"
  )

# make new column for focal or anc_proxy, if species is ilMecLysi212 and events ending in Q1 -> focal
#if species is ilMecLysi212 and events ending in Q2 -> anc_proxy and vice versa for ilMecPoly1, if starts with shared then shared

df <- df %>%
  mutate(
    event_category = case_when(
      species == "ilMecPoly1" & endsWith(type, "Q1") ~ "focal",
      species == "ilMecPoly1" & endsWith(type, "Q2") ~ "anc_proxy",
      species == "ilMecLysi212" & endsWith(type, "Q2") ~ "focal",
      species == "ilMecLysi212" & endsWith(type, "Q1") ~ "anc_proxy",
      grepl("shared", type) ~ "shared_event",
      TRUE ~ "other"
    )
  )


head(df)
str(df)
# merge all events that have the same chr, start, end, type, event_parent_node and species, and if they occur in this genome or the anc proxy.
# merge brp id in one column and if the type is the same or different, if different then list all types involved
df_2 <- df %>%
  group_by(chr, start, end, species, event_parent_node, event_category) %>%
  summarise(
    all_brp_ids = paste(all_brp_ids, collapse = ","),
    #if type is fission_Q1 and fusion_Q1
    type = case_when(
      setequal(sort(unique(type)), c("fission_Q1", "fusion_Q1")) ~ "fission-fusion_Q1",
      setequal(sort(unique(type)), c("fission_Q2", "fusion_Q2")) ~ "fission-fusion_Q2",
      n_distinct(type) == 1 ~ first(type),
      TRUE ~ paste(sort(unique(type)), collapse = ",")
    ),
    .groups = "drop"
  ) %>%
  mutate(type = sub("(_Q1|_Q2)$", "", type))

head(df_2)
str(df_2)

unique(df_2$event_category)

df_2 %>% filter(event_category != "anc_proxy") %>% group_by(species, type, event_parent_node, event_category) %>% summarise(n = n())

sink("tables/anc_breakpoint_counts_by_species_and_type_notes.txt", append = TRUE)
print("Summary of breakpoint counts after determining most recent node and all nodes involved:")
df_2 %>% group_by(species, type, event_parent_node) %>%summarise(n = n())
sink()


write.table(df_2, "tables/anc_breakpoints_formatted.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE) 

write.table(df_2 %>% filter(event_category != "anc_proxy") %>% group_by(species, type, event_parent_node, event_category) %>% summarise(n = n()),
            "tables/anc_breakpoint_counts_by_species_and_type.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

df_2 %>% filter(event_category != "anc_proxy") %>% 
  ggplot() +
  geom_bar(aes(x = event_parent_node, fill = type), position = "dodge") +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Counts of breakpoint events by species, type, and parent node",
       x = "Parent node",
       y = "Count of events",
       fill = "Event type")


ggsave("plots/anc_breakpoint_counts_by_species_and_type_before_correction.pdf", width = 8, height = 6)


###

#remove lines where event_parent_node is N1 and type is fission, 

# as these are likely misclassified, all older event should be shared unless recurrent events
# after checking synteny with close relatives these are likely misscalssifications and will be removed

df_2 %>%
  filter((event_parent_node %in% c("N1") & (event_category == "focal" | event_category == "anc_proxy")))


df_2 <- df_2 %>%
  mutate(
    type = case_when(
      event_parent_node == "N1" & (event_category == "focal" | event_category == "anc_proxy") ~ "Complex",
      TRUE ~ type
    )
  )

View(df_2)
# Missasigned shared fusion at N3, likely a misclassification, 
# uncertain the fusion may be shared with an inversion or two separate events
df_2 %>% filter((event_parent_node == "N3" & type == "shared_fusion"))

df_2 <- df_2 %>%
  mutate(
    type = case_when(
      event_parent_node == "N3" & type == "shared_fusion" ~ "Complex",
      TRUE ~ type
    )
  )


write.table(df_2, "tables/anc_breakpoints_formatted_after_correction.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE) 

write.table(df_2 %>% filter(event_category != "anc_proxy") %>% 
              group_by(species, type, event_parent_node, event_category) %>% 
              summarise(n = n()),
            "tables/anc_breakpoint_counts_by_species_and_type_after_correction.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


df_2 %>% filter(event_category != "anc_proxy") %>% 
  ggplot() +
  geom_bar(aes(x = event_parent_node, fill = type), position = "dodge") +
  facet_wrap(~species) +
  theme_bw() +
  labs(title = "Counts of breakpoint events by species, type, and parent node",
       x = "Parent node",
       y = "Count of events",
       fill = "Event type")


ggsave("plots/anc_breakpoint_counts_by_species_and_type_after_correction.pdf", width = 8, height = 6)

