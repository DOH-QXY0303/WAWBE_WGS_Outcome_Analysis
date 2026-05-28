###############################################
######## 1) LOAD REQUIRED LIBRARIES ###########
###############################################

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(DBI)
library(odbc)
library(keyring)
library(purrr)
library(writexl)

#Connect to LIMS
con <- dbConnect(
  odbc(),
  Driver   = "ODBC Driver 18 for SQL Server",
  Server   = "DOH01DBTUMP22,9799",
  Database = "LIMS_DATA",
  TrustServerCertificate = "Yes",
  Trusted_Connection = "Yes"
  )

lims_data <- dbGetQuery(
  con,
  "SELECT * FROM [vz_Epi_ELS_Multiplex_SARS-Flu-RSV-Mpox_dPCR]"
)

###############################################
### 2)IMPORT EXCEL, APPEND TABS, CLEAN DATA ###
###############################################

# Step 1: Open file browser and save the selected Excel file path as a character string
file_path <- file.choose()

# Step 2: Get all worksheet/tab names from the selected Excel file
sheets <- excel_sheets("/mnt/c/Users/qxy0303/scratch/WAWBE_outcome_analysis/Run Stats.xlsx")

# Step 3: Loop through each sheet:
#   - Read the sheet into a dataframe
#   - Add a new column (IC_Group) containing the sheet name
#   - Store each sheet dataframe in a list
# Step 4: Combine all sheet dataframes into one large dataframe
df <- lapply(sheets, function(sheet_name) {
  read_excel(file_path, sheet = sheet_name) %>%
    mutate(IC_Group = sheet_name)
}) %>%
  bind_rows()


# Step 5: Change column names 
colnames(df)[colnames(df) == "Key_ID"] <- "PHLAccessionNumber"
colnames(df)[colnames(df) == ">10X Coverage"] <- "Coverage_10X"

# Step 6: Merge WWTP names from LIMS to df by PHAccessionNumber
df = merge(df, lims_data[,c("PHLAccessionNumber","WWTPName")], by="PHLAccessionNumber")


# Step 7: Remove columns
df <- df %>%
  select(-all_of(c(
    "0x Coverage",
    "...10",
    "...11",
    "#A",
    "#B",
    "#A/B",
    "#C",
    "% Failed",
    "% PASS (A)",
    "% PASS (A/B)",
    "Column1",
    "Library Prep Input Copies"
  )))

# Step 8: Clean data frame: remove rows where all columns are NA and where IC_Group is "Val Table"
df <- df %>%
  filter(
    rowSums(is.na(.)) != ncol(.),   # remove rows that are entirely NA
    IC_Group != "Val Table"         # remove unwanted group
  )

#TESTING THE CORRELATION BETWEEN IC AND >10X COVERAGE#

# Correlation coefficient only
cor(df$`Input Copies`, df$Coverage_10X, method = "pearson", use = "complete.obs")

# Correlation coefficient with significance test (p-value)
cor.test(df$`Input Copies`, df$Coverage_10X, method = "pearson")

###############################################
##### MAKE SCATTER PLOTS BOX PLOTS GRAPHS #####
###############################################

#Correlation Scatter plot
ggplot(df, aes(x = `Input Copies`, y = Coverage_10X)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  geom_hline(yintercept = 90, color = "red", linewidth = 1.5) +
  labs(
    title = "Correlation Between Concentration and >10X Coverage",
    x = "Number of Input Copies",
    y = "Samples with >10X Coverage"
  ) +
  scale_x_continuous(
    limits = c(0, 200)
  ) +
  scale_y_continuous(
    breaks = c (0, 25, 50, 75, 90, 100))
                     

#Boxplot of Input Copies
ggplot(df, aes(x = WWTPName, y = `Input Copies`, fill = WWTPName)) +
  geom_boxplot() +
  ggtitle("Input Copies by WWTP") +
  scale_y_continuous(limits = c(0, 200)) +
  scale_x_discrete(
    labels = function(x)
      gsub(
        "WWTP|Water|Reclamation|Facility|Treatment|Plant|Wastewater|WRP|Clean|\\(SCTP\\)|Influent|Regional|and|Sewage",
        "",
        x
      )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )

#Boxplot of >10X Coverage
ggplot(df, aes(x = WWTPName, y = Coverage_10X, fill = WWTPName)) +
  geom_boxplot() +
  ggtitle("Samples with >10X by WWTP") +
  scale_x_discrete(
    labels = function(x)
      gsub(
        "WWTP|Water|Reclamation|Facility|Treatment|Plant|Wastewater|WRP|Clean|\\(SCTP\\)|Influent|Regional|and|Sewage",
        "",
        x
      )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )

#Table
quality_summary <- df %>%
  mutate(
    Quality = case_when(
      Coverage_10X >= 90 ~ "A",
      Coverage_10X >= 60 & Coverage_10X < 90 ~ "B",
      Coverage_10X < 60 ~ "Fail",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(WWTPName) %>%
  summarise(
    Total_Samples = n(),
    
    Percent_A = round(mean(Quality == "A", na.rm = TRUE) * 100, 1),
    
    Percent_B = round(mean(Quality == "B", na.rm = TRUE) * 100, 1),
    
    Percent_Fail = round(mean(Quality == "Fail", na.rm = TRUE) * 100, 1),
    
    Avg_Input_Copies = round(mean(`Input Copies`, na.rm = TRUE), 2)
  ) %>%
  arrange(desc(Percent_A))

qs <- quality_summary

write.csv(qs, file = "WWTP by Sample Quality", row.names = FALSE)
write_xlsx(qs, "quality_summary.xlsx")

###############################################
############# 7) CREATE TABLE #################
###############################################

#creates groupings for input copies and removes QCD with N/As

df2 <- df %>%
  filter(!is.na(QCD)) %>%
  mutate(
    IC_Group = case_when(
      `Input Copies` >= 0   & `Input Copies` <= 4.9    ~ "IC 0-4",
      `Input Copies` >= 5   & `Input Copies` <= 9.9    ~ "IC 5-9",
      `Input Copies` >= 10  & `Input Copies` <= 19.9   ~ "IC 10-19",
      `Input Copies` >= 20  & `Input Copies` <= 29.9   ~ "IC 20-29",
      `Input Copies` >= 30  & `Input Copies` <= 39.9   ~ "IC 30-39",
      `Input Copies` >= 40  & `Input Copies` <= 49.9   ~ "IC 40-49",
      `Input Copies` >= 50  & `Input Copies` <= 100.9  ~ "IC 50-100",
      `Input Copies` >= 101 & `Input Copies` <= 1000 ~ "IC 101-1000",
      TRUE ~ "Other"
    )
  )

#creates pivot table?
Ic_Group_Table <- df2 %>%
  
  group_by(WWTPName, IC_Group, QCD) %>%
  summarise(n = n(), .groups = "drop") %>%
  
  pivot_wider(
    names_from = c(IC_Group, QCD),
    values_from = n,
    values_fill = 0
  )

write.csv(df2, "WWTP_QCD_Data_2.csv", row.names = FALSE)



