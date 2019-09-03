
clean_data(file.path(getwd(), "/experiments/Jon_Uq2_unc45/1239Shp15a_NP/PID summary"))

clean_data(file.path(getwd(), "/experiments/Jon_Uq2_unc45/JK11_Uq2_78/PID summary"))

working_path <- "experiments/Jon_Uq2_unc45/"

data <- bind_rows(read_csv(file.path(getwd(), working_path, "JK11_Uq2_78/cleaned_CDR3s.csv")),
                  read_csv(file.path(getwd(), working_path, "1239Shp15a_NP/cleaned_CDR3s.csv")))

data_clustered <- data %>%
  left_join(all_clusters)
