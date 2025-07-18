# Quick test to check oracle output structure
library("dplyr")
library("hubUtils")
library("hubData")

hub_path <- here::here("sismid-ili-forecasting-sandbox")

oracle_output <- connect_target_oracle_output(hub_path) |> 
  collect()

print("Oracle output structure:")
print(names(oracle_output))
print(head(oracle_output))