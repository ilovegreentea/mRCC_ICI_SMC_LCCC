library(maftools)
library(dplyr)
library(readr)

sample_annot <- tribble(
  ~Sample_ID, ~Patient_ID, ~Timepoint, ~Response,
  "D_23_00810_TS_WES", "RCC5", "T2", "Responder",
  "D_23_00811_TS_WES", "RCC7", "T0", "Responder",
  "D_23_00813_TS_WES", "RCC8", "T2", "Responder",
  "D_24_00121_TS_WES", "RCC10", "T2", "Non_responder",
  "D_23_00396_TS_WES", "RCC12", "T0", "Non_responder",
  "D_23_00397_TS_WES", "RCC13", "T0", "Non_responder",
  "D_23_00399_TS_WES", "RCC14", "T2", "Non_responder",
  "D_23_00401_TS_WES", "RCC15", "T0", "Non_responder",
  "D_23_00402_TS_WES", "RCC16", "T0", "Responder",
  "D_23_00403_TS_WES", "RCC17", "T0", "Responder",
  "D_23_00404_TS_WES", "RCC18", "T0", "Responder",
  "D_23_00405_TS_WES", "RCC19", "T0", "Responder",
  "D_23_00407_TS_WES", "RCC20", "T2", "Non_responder",
  "D_23_00408_TS_WES", "RCC21", "T0", "Non_responder",
  "D_23_00814_TS_WES", "RCC22", "T0", "Responder",
  "D_23_00815_TS_WES", "RCC23", "T2", "Responder",
  "D_23_00816_TS_WES", "RCC24", "T0", "Non_responder",
  "D_23_00817_TS_WES", "RCC25", "T2", "Non_responder",
  "D_23_00818_TS_WES", "RCC26", "T0", "Responder"
)
data_dir <- "C:/Users/gkmin/OneDrive/desktop/main_project/IO/script/hg38_maf/"
maf_files <- list.files(data_dir, full.names = TRUE, pattern = "*.maf$")

maf_list <- lapply(maf_files, read.maf)

merged_maf <- merge_mafs(maf_list)
library(data.table)

sample_annot_dt <- as.data.table(sample_annot)

setnames(sample_annot_dt, "Sample_ID", "Tumor_Sample_Barcode")

merged_maf@clinical.data <- sample_annot_dt
maf_filtered <- merged_maf

filtered_data <- maf_filtered@data %>%
  filter(
    n_alt_count == 0,            
    t_depth >= 20,               
    t_alt_count / t_depth >= 0.2, 
    FILTER == "PASS"            
  )

filtered_data <- filtered_data %>%
  filter(Variant_Classification %in% c(
    "Missense_Mutation",
    "Splice_Site",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Ins",     
    "In_Frame_Del",    
    "Nonsense_Mutation"
  ))

maf_filtered@data <- filtered_data

oncoplot(maf_filtered, clinicalFeatures = "Response", sortByAnnotation = TRUE)

library(maftools)
library(dplyr)

setdiff(original_samples, remaining_samples)

filtered_df <- as.data.frame(
  maf_filtered@data %>%
    filter(!grepl("^MT-", Hugo_Symbol))
)


maf_filtered_nomt <- read.maf(
  maf = filtered_df,
  clinicalData = maf_filtered@clinical.data
)
maf_filtered_nomt@clinical.data

patient_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(unique(maf_filtered_nomt@clinical.data$Patient_ID))),
  unique(maf_filtered_nomt@clinical.data$Patient_ID)
)
maf_filtered_nomt@clinical.data
oncoplot(
  maf = maf_filtered_nomt,
  clinicalFeatures = c("Response", "Timepoint", "Patient_ID"),
  sortByAnnotation = TRUE,
  annotationColor = list(
    Response = c(
      "Responder" = "#00BFC4",
      "Non_responder" = "#F8766D"
    ),
    Timepoint = c(
      "T0" = "#7CAE00",
      "T2" = "#C77CFF"
    ),
    Patient_ID = patient_colors
  )
)
