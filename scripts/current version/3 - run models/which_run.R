library(tidyverse)
library(here)

# results_dir <- here("loon_results")

# Step 1: List files matching the more general pattern
files <- list.files(here('results'), pattern = "^(highout|lowout|medout)-\\d+-\\d+-\\d+\\.RDS$", full.names = FALSE)

# Step 2: Extract the numbers as before
extract_numbers <- function(filename) {
  # Remove ".RDS"
  name <- sub("\\.RDS$", "", filename)
  # Split by "-"
  parts <- strsplit(name, "-")[[1]]
  # Get last three elements
  nums <- tail(parts, 3)
  # Convert to integer
  as.integer(nums)
}

# Step 3: Build the dataframe
numbers <- t(sapply(files, extract_numbers))
df <- data.frame(
  file = files,
  first = numbers[,1],
  middle = numbers[,2],
  last = numbers[,3]
)
rownames(df) <- NULL

print(df)

has_run <- df %>% 
  mutate(
    traj = case_when(
      str_detect(file, "high") ~ "high",
      str_detect(file, "low") ~ "low",
      TRUE ~ "med"
    ), 
    .before = 1
  ) %>% 
  select(-file)
were_queued <- expand_grid(
  traj = c("high", "med", "low"),
  first = 1:3, 
  middle = 1:12,
  last = 1:48
)
hasnt_run <- anti_join(were_queued, has_run) %>% 
  # because code runs over all trajectories and data scenarios
  select(-c(traj, last)) %>% 
  distinct()
new_sims_per <- hasnt_run$middle %>% unique()
new_scenarios_picked <- hasnt_run$first %>% unique()
