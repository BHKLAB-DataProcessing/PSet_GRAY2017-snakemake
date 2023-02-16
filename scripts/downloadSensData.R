library("data.table")
library("matrixStats")
library(abind)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- paste0(args[1], "data")
download_dir <- paste0(args[1], "download")
processed_dir <- paste0(args[1], "processed")

# data_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2017/data"
# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2017/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/getGRAY2017/processed"

concentrations.no <- 9

options(stringsAsFactors = FALSE)
gray.raw.drug.sensitivity <- fread(file = file.path(data_dir, "Gray_data_raw_dose_response.csv"))
gray.conc <- fread(file.path(data_dir, "Gray_drug_conc.csv"))

gray.raw.drug.sensitivity[
  ,
  background := rowMedians(data.matrix(gray.raw.drug.sensitivity[, grep("^background_od",
    colnames(gray.raw.drug.sensitivity),
    v = TRUE
  ), with = FALSE]))
]

gray.raw.drug.sensitivity[
  ,
  ctrl := rowMedians(data.matrix(gray.raw.drug.sensitivity[, grep("^od0\\.",
    colnames(gray.raw.drug.sensitivity),
    v = TRUE
  ), with = FALSE])) - background
]


for (i in seq_len(concentrations.no)) {
  curCol <- paste0("response", i)

  curDoseCol <- grep(paste0("^od", i, "\\."), colnames(gray.raw.drug.sensitivity))
  gray.raw.drug.sensitivity[
    ,
    (curCol) := (rowMedians(data.matrix(gray.raw.drug.sensitivity[, curDoseCol,
      with = FALSE
    ])) - background) / ctrl * 100
  ]
}

gray.raw.drug.sensitivity <- gray.raw.drug.sensitivity[, -grep("^od",
  colnames(gray.raw.drug.sensitivity),
  v = TRUE
), with = FALSE]

# gray.raw.drug.sensitivity <- gray.raw.drug.sensitivity[,-grep("*background*",
#                                                               colnames(gray.raw.drug.sensitivity),
#                                                               v=TRUE), with=FALSE]
# gray.raw.drug.sensitivity <- gray.raw.drug.sensitivity[,-c("T0_median_od")]



gray.raw.drug.sensitivity <- merge(
  x = gray.raw.drug.sensitivity,
  y = gray.conc, by = c("drug_group_id", "drug"), all.x = TRUE
)



gray.raw.drug.sensitivity[, exp_id := paste(drug, cellline, ID.x, sep = "_")]

conc.cols <- grep("^c[0-9]+", colnames(gray.raw.drug.sensitivity), v = TRUE)

gray.raw.drug.sensitivity[, (conc.cols) := lapply(.SD, function(x) {
  return(x * 1e6)
}), .SDcols = conc.cols]

gray.raw.drug.sensitivity[drug == "2C4", (conc.cols) := .SD / 388.8, .SDcols = conc.cols] ## molecular weight from Pubchem, retrieved Jul 8, 2020
gray.raw.drug.sensitivity[drug == "2C4", units := "M"] # removing one drug in mg/ml for now.



# gray.conc.m <- data.table::melt(gray.raw.drug.sensitivity[,c("exp_id", conc.cols), with=FALSE], id.vars="exp_id")

# colnames(gray.conc.m) <- c("exp_id", "DoseNum", "DoseValue")

# gray.conc.m[,DoseNum := gsub(DoseNum, pat="c", rep="")]

response.cols <- grep("response", x = colnames(gray.raw.drug.sensitivity), v = TRUE)

raw.sensitivity <- abind(data.matrix(gray.raw.drug.sensitivity[, ..conc.cols]), data.matrix(gray.raw.drug.sensitivity[, ..response.cols]), along = 3)

rownames(raw.sensitivity) <- gray.raw.drug.sensitivity[["exp_id"]]

colnames(raw.sensitivity) <- paste0("dose", seq_len(concentrations.no))

dimnames(raw.sensitivity)[[3]] <- c("Dose", "Viability")

sensitivity.info <- as.data.frame(gray.raw.drug.sensitivity[, .(drug, cellline, drug_plate_id, T0_plate_id, T0_median_od, T0_background_od1, T0_background_od2, units)])
colnames(sensitivity.info)[1:2] <- c("drugid", "cellid")
sensitivity.info$nbr.conc.tested <- 9
sensitivity.info$min.Dose.uM <- apply(raw.sensitivity[, , "Dose"], 1, min)
sensitivity.info$max.Dose.uM <- apply(raw.sensitivity[, , "Dose"], 1, max)
rownames(sensitivity.info) <- rownames(raw.sensitivity)
save(raw.sensitivity, sensitivity.info, file = file.path(processed_dir, "drug_norm_post.RData"))

## create sensitivity slices


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity) / 1000))

dir.create(file.path(processed_dir, "slices"))

for (i in seq_along(raw.sensitivity.x)) {
  slce <- raw.sensitivity[raw.sensitivity.x[[i]], , ]
  saveRDS(slce, file = paste0(processed_dir, "/slices/gray2017_raw_sens_", i, ".rds"))
}
zip(zipfile=file.path(processed_dir, 'raw_sense_slices.zip'), files=list.files(file.path(processed_dir, 'slices'), full.names = TRUE))
unlink(file.path(processed_dir, 'slices'), recursive=TRUE)
