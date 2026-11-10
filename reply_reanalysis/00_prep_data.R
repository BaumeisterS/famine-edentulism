# =====================================================================
#  00_prep_data.R  —  Re-analysis of Baumeister et al. 2026, J Clin
#  Periodontol 53:1673-80 (doi:10.1111/jcpe.70190), prepared for the
#  authors' reply to correspondence on the article.
#
#  Builds slim, local analysis files (individual-level; NEVER share or
#  commit — they stay in data_local/):
#   * clhls_slim.rds  : CLHLS 2008 baseline with the RAW natural-tooth
#                       count (g21), birth province (a41), birth year, age,
#                       sex -> lets us re-code the edentulism threshold.
#   * charls_did.rds  : CHARLS wave-2 edentulism + birth province
#                       (Life-History Residence module, first spell),
#                       exactly as in the published analysis code.
#  Variable definitions copy the published script
#  (github/famine_edentulism_analysis.R) line for line.
# =====================================================================
suppressPackageStartupMessages({ library(haven); library(dplyr) })
DAT <- Sys.getenv("FAMINE_DATA_DIR", unset = "c:/arbeit/Greifswald/Daten/")
LOC <- "data_local/"
t0 <- Sys.time()

# ---- CLHLS (Henan famine 1942-43) -----------------------------------
cl <- read_dta(paste0(DAT, "CLHLS/clhls_2008_2018_longitudinal_dataset_released_version1.dta"),
               col_select = any_of(c("id","v_bthyr","trueage","a1","a41","g21","f66","f85","f86")))
clhls <- cl |>
  transmute(byear = as.numeric(v_bthyr),
            age   = as.numeric(trueage),
            male  = case_when(a1 == 1 ~ 1L, a1 == 2 ~ 0L),
            g21   = as.numeric(g21),
            teeth = if_else(g21 > 32, NA_real_, g21),          # natural teeth 0-32
            region = as.character(as_factor(a41)),
            father_occupation = case_when(f85 %in% c(0,1,3) ~ 0L, TRUE ~ 1L),
            father_education  = case_when((f86 >= 0 & f86 < 10) | f86 %in% c(88,99) ~ 1L,
                                           f86 >  9 & f86 < 12 ~ 2L, f86 > 11 ~ 3L))
saveRDS(clhls, paste0(LOC, "clhls_slim.rds"))
cat("CLHLS rows:", nrow(clhls), "| teeth non-missing:", sum(!is.na(clhls$teeth)), "\n")
cat("CLHLS natural-teeth distribution (0..6 and >6):\n")
print(table(cut(clhls$teeth, c(-1,0,1,2,3,4,5,6,32), labels=c("0","1","2","3","4","5","6","7-32")), useNA="ifany"))

# ---- CHARLS (Great Famine 1959-61) ----------------------------------
gb <- c('11'='Beijing','12'='Tianjin','13'='Hebei','14'='Shanxi','15'='Inner Mongolia',
        '21'='Liaoning','22'='Jilin','23'='Heilongjiang','31'='Shanghai','32'='Jiangsu',
        '33'='Zhejiang','34'='Anhui','35'='Fujian','36'='Jiangxi','37'='Shandong',
        '41'='Henan','42'='Hubei','43'='Hunan','44'='Guangdong','45'='Guangxi',
        '46'='Hainan','50'='Chongqing','51'='Sichuan','52'='Guizhou','53'='Yunnan',
        '54'='Tibet','61'='Shaanxi','62'='Gansu','63'='Qinghai','64'='Ningxia','65'='Xinjiang')
res <- read_dta(paste0(DAT, "CHARLS/2014_life_history/Residence.dta"),
                col_select = any_of(c("ID","t002_psu_1_")))
res$psu1 <- as.numeric(res$t002_psu_1_)
bp <- tibble(ID = res$ID,
             birthprov = unname(gb[as.character(ifelse(res$psu1 >= 1100, floor(res$psu1/100), NA))]))
h  <- read_dta(paste0(DAT, "CHARLS/harmonized/H_CHARLS_D_Data.dta"),
               col_select = any_of(c("ID","rabyear","r1agey","ragender")))
ed <- read_dta(paste0(DAT, "CHARLS/w2/Health_Status_and_Functioning.dta"),
               col_select = any_of(c("ID","da040")))
charls <- h |>
  transmute(ID, byear = as.numeric(rabyear), age = as.numeric(r1agey),
            male = if_else(ragender == 1, 1L, 0L)) |>
  left_join(ed |> transmute(ID, edent = case_when(da040 == 1 ~ 1L, da040 == 2 ~ 0L)), "ID") |>
  left_join(bp, "ID")
saveRDS(charls, paste0(LOC, "charls_did.rds"))
cat("\nCHARLS rows:", nrow(charls), "| with birth province:", sum(!is.na(charls$birthprov)), "\n")

# cross-check against the birth-province file used during the revision
old <- readRDS(paste0(LOC, "charls_birthprov.rds"))
chk <- inner_join(bp, old |> mutate(birthprov_old = unname(birthprov)) |> select(ID, birthprov_old), "ID")
cat("Birth-province agreement with revision file:",
    round(100*mean(chk$birthprov == chk$birthprov_old, na.rm = TRUE), 2), "% of",
    sum(!is.na(chk$birthprov) & !is.na(chk$birthprov_old)), "matched IDs\n")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s\n")
