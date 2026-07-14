# =====================================================================
#  Childhood famine exposure and late-life edentulism:
#  evidence from natural experiments
#
#  Manuscript 9424511, Journal of Clinical Periodontology (revised)
#  Baumeister S-E, McGrath C, Kocher T, Nascimento GG
#
#  This single script reproduces every number, table, and figure in the
#  manuscript and supplement from harmonised individual-level data on
#  ten historical famines across Asia, Europe and the Americas.
#
#  DESIGN RATIONALE
#  ----------------
#  Birth-year-based famine exposure confounds exposure with AGE AT
#  SURVEY (the famine-exposed cohort is older, and edentulism rises
#  steeply with age). The crude reduced-form RR collapses from 1.69
#  to 1.07 after age adjustment. We therefore use quasi-experimental
#  designs that separate famine exposure from the age gradient:
#
#   PRIMARY  — Regression Discontinuity in Time (RDiT; Hausman &
#              Rapson 2018) with bias-corrected CIs (rdrobust;
#              Calonico et al. 2014), honest CIs for the discrete
#              running variable (RDHonest; Armstrong & Kolesár 2018),
#              and Donut RDiT. Applied to all ten cohorts.
#   SECONDARY— Within-country region-of-birth × famine-intensity DiD
#              (Chen & Zhou 2007) in the two Chinese cohorts (CHARLS
#              Great Famine 1959-61; CLHLS Henan famine 1942-43).
#              Unadjusted TWFE primary; adjusted + BJS imputation as
#              sensitivity. Single-event 2×2 design → TWFE unbiased
#              (Baker et al. 2026).
#   SECONDARY— Interrupted time-series (segmented regression with
#              level shift and slope change at the famine cohort).
#   VALIDITY — Three negative controls (placebo cohort, placebo
#              treatment window age 11-15, negative-control outcome).
#
#  SENSITIVITY (SUPPLEMENT)
#   S1: Descriptive statistics by study and cohort
#   S2: Reduced-form RR before and after age adjustment (1.69→1.07)
#   S3: Leave-one-out (incl. mortality-selected Mexico)
#   S4: Famine-duration stratification (≤2 vs >2 years)
#
#  Childhood exposure is located by PROVINCE OF BIRTH (not survey
#  province). CHARLS: Life-History Residence module t002_psu_1_.
#  CLHLS: a41 "province of birth".
#
#  Risk differences (percentage points, linear probability models) are
#  reported throughout; the reduced-form is also reported as a Poisson
#  risk ratio for comparability with the original submission.
# =====================================================================

# ---- setup -------------------------------------------------------
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr)
  library(fixest); library(rdrobust); library(rddensity); library(RDHonest)
  library(didimputation)
  library(sandwich); library(lmtest); library(broom)
  library(ggplot2); library(forcats)
})

DAT <- Sys.getenv("FAMINE_DATA_DIR",
                   unset = "c:/arbeit/Greifswald/Daten/")
OUT <- "."                                      # tables/ and figures/ live next to this script
FIG <- paste0(OUT, "/figures/"); TAB <- paste0(OUT, "/tables/")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base = 10)
  theme_classic(base_size = base) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

study_labels <- c("1"="CHARLS","2"="CLHLS","3"="LASI","4"="SAGE-India",
                  "5"="SAGE-Mexico","6"="SAGE-Russia","7"="SHARE-Netherlands",
                  "8"="SHARE-Germany","9"="SHARE-Austria","10"="SHARE-Spain")

# =====================================================================
# SECTION 1 — DATA PREPARATION  (10 studies)
# =====================================================================
# Cached harmonised file speeds re-runs. Delete .rds to force rebuild.
.cache <- paste0(OUT, "/famine_edent_pooled.rds")
if (file.exists(.cache)) {
  cat("== Cached pooled data found — loading ==\n")
  pooled <- readRDS(.cache)
} else {

# -- 1.1 CHARLS (Chinese Great Famine 1959-61) ----
charls <- read_dta(paste0(DAT, "CHARLS/harmonized/H_CHARLS_D_Data.dta")) |>
  mutate(famine_cohort_1_5_itt = case_when(rabyear %in% 1953:1960 ~ 1L,
                                            rabyear %in% 1961:1969 ~ 0L),
         age  = r1agey,
         male = if_else(ragender == 1, 1L, 0L, missing = NA_integer_))

lh_charls <- read_dta(paste0(DAT, "CHARLS/2014_life_history/H_CHARLS_LH_a.dta")) |>
  mutate(father_occupation = case_when(radadoccup_c == 1 ~ 1L, radadoccup_c == 2 ~ 0L)) |>
  select(ID, father_occupation)

cs_charls <- read_dta(paste0(DAT, "CHARLS/2014_life_history/Family_Information.dta")) |>
  mutate(childhood_starvation = case_when(c3_a == 1 ~ 1L, c3_a == 2 ~ 0L)) |>
  select(ID, childhood_starvation)

edent_charls <- read_dta(paste0(DAT, "CHARLS/w2/Health_Status_and_Functioning.dta")) |>
  mutate(edentulism_2 = case_when(da040 == 1 ~ 1L, da040 == 2 ~ 0L)) |>
  select(ID, edentulism_2)

charls_harm <- charls |>
  left_join(lh_charls,    by = "ID") |>
  left_join(cs_charls,    by = "ID") |>
  left_join(edent_charls, by = "ID") |>
  transmute(study = 1L, country = 1L, continent = 1L,
            edent = edentulism_2, treatment = famine_cohort_1_5_itt,
            byear = rabyear, famine_start = 1959L, famine_end = 1961L,
            age, male, father_occupation, childhood_starvation)

# -- 1.2 CLHLS (Henan famine 1942-43) ----
clhls_harm <- read_dta(paste0(DAT,
  "CLHLS/clhls_2008_2018_longitudinal_dataset_released_version1.dta")) |>
  mutate(childhood_starvation = case_when(f66 == 1 ~ 1L, f66 == 2 ~ 0L),
         famine_cohort_1_5_itt = case_when(v_bthyr %in% 1937:1942 ~ 1L,
                                            v_bthyr %in% 1943:1948 ~ 0L),
         num_teeth = if_else(g21 > 32, NA_real_, g21),
         num_missingteeth28 = case_when((32 - num_teeth) %in% 28:32 ~ 1L,
                                         (32 - num_teeth) %in%  0:27 ~ 0L),
         male = case_when(a1 == 1 ~ 1L, a1 == 2 ~ 0L),
         age  = trueage,
         father_occupation = case_when(f85 %in% c(0,1,3) ~ 0L, TRUE ~ 1L),
         father_education  = case_when((f86 >= 0 & f86 < 10) | f86 %in% c(88,99) ~ 1L,
                                        f86 >  9 & f86 < 12                       ~ 2L,
                                        f86 > 11                                  ~ 3L)) |>
  transmute(study = 2L, country = 1L, continent = 1L,
            edent = num_missingteeth28, treatment = famine_cohort_1_5_itt,
            byear = v_bthyr, famine_start = 1942L, famine_end = 1945L,
            age, male, father_occupation, father_education, childhood_starvation)

# -- 1.3 LASI (Bengal famine 1943-44) ----
lasi_harm <- read_dta(paste0(DAT, "LASI/H_LASI_a3.dta")) |>
  transmute(study = 3L, country = 2L, continent = 1L,
            edent = r1noteeth,
            treatment = case_when(rabyear %in% 1938:1943 ~ 1L,
                                  rabyear %in% 1944:1945 ~ 0L),
            byear = rabyear, famine_start = 1943L, famine_end = 1944L,
            age = r1agey,
            male = if_else(ragender == 1, 1L, 0L, missing = NA_integer_),
            father_education = radadeducl)

# -- 1.4-1.6 SAGE (India, Mexico, Russia) ----
sage_raw <- read_dta(paste0(DAT, "SAGE/wave0/ICPSR_28502/DS0001/28502-0001-Data.dta"))

make_sage <- function(cntry, fam_born, post_born, fam_start, fam_end,
                      study_id, country_id, cont_id) {
  sage_raw |> filter(COUNTRY == cntry) |>
    transmute(study = study_id, country = country_id, continent = cont_id,
              edent = case_when(Q6757 == 1 ~ 1L, Q6757 == 5 ~ 0L),
              age   = Q1002, byear = 2003L - as.integer(Q1002),
              treatment = case_when(byear %in% fam_born  ~ 1L,
                                    byear %in% post_born ~ 0L),
              male = case_when(Q1001 == 2 ~ 1L, Q1001 == 1 ~ 0L),
              famine_start = fam_start, famine_end = fam_end)
}
sage_india_harm  <- make_sage("IND", 1938:1943, 1944:1948, 1943L, 1944L, 4L, 2L, 1L)
sage_mexico_harm <- make_sage("MEX", 1909:1916, 1917:1920, 1914L, 1917L, 5L, 3L, 2L)
sage_russia_harm <- make_sage("RUS", 1940:1946, 1947:1951, 1946L, 1947L, 6L, 4L, 1L)

# -- 1.7-1.10 SHARE (NL, DE, AT, ES) ----
share_raw <- read_dta(paste0(DAT,
  "SHARE/GateWay_Harmonized_SHARE_Version9/GH_SHARE_g.dta"))

make_share <- function(ccode, fam_born, post_born, fam_start, fam_end,
                       study_id, country_id, cont_id) {
  share_raw |> filter(country == ccode) |>
    transmute(study = study_id, country = country_id, continent = cont_id,
              edent = r5noteeth,
              treatment = case_when(rabyear %in% fam_born  ~ 1L,
                                    rabyear %in% post_born ~ 0L),
              byear = rabyear, famine_start = fam_start, famine_end = fam_end,
              age  = coalesce(r5agey, r4agey, r2agey, r1agey),
              male = if_else(ragender == 1, 1L, 0L, missing = NA_integer_),
              father_education = radadeducl)
}
share_dutch_harm   <- make_share(14, 1939:1944, 1946:1950, 1944L, 1945L,  7L, 6L, 3L)
share_germany_harm <- make_share(12, 1940:1947, 1948:1952, 1945L, 1948L,  8L, 7L, 3L)
share_austria_harm <- make_share(11, 1940:1947, 1948:1952, 1945L, 1948L,  9L, 5L, 3L)
share_spain_harm   <- make_share(15, 1934:1941, 1942:1948, 1939L, 1942L, 10L, 8L, 3L)

# -- Pool and derive analysis variables ----
pooled <- bind_rows(
  charls_harm, clhls_harm, lasi_harm,
  sage_india_harm, sage_mexico_harm, sage_russia_harm,
  share_dutch_harm, share_germany_harm, share_austria_harm,
  share_spain_harm) |>
  mutate(
    study_f    = factor(study, levels = 1:10, labels = study_labels),
    famine_dur = famine_end - famine_start + 1,
    short_famine = as.integer(famine_dur <= 2),
    event_time = byear - famine_start,
    cohort3 = case_when(
      byear <  famine_start - 5                          ~ 0L,
      byear >= famine_start - 5 & byear <= famine_end - 1 ~ 1L,
      byear >  famine_end                                 ~ 2L),
    cohort3_f = factor(cohort3, 0:2,
                       c("Pre-famine","Famine-exposed","Post-famine")),
    father_ses = case_when(
      !is.na(father_education) & father_education == 1 ~ 1L,
      !is.na(father_education) & father_education == 2 ~ 2L,
      !is.na(father_education) & father_education == 3 ~ 3L,
      !is.na(father_occupation) & father_occupation == 1 ~ 1L,
      !is.na(father_occupation) & father_occupation == 0 ~ 2L,
      TRUE ~ 1L),
    father_ses_f = factor(father_ses, 1:3, c("low/missing","medium","high")),
    age_c  = age - 70, age_c2 = age_c^2, age_c3 = age_c^3,
    male   = as.integer(male),
    edent  = as.integer(edent))

cat("\n== POOLED N =", nrow(pooled), "==\n")
print(with(pooled, table(study_f, cohort3_f, useNA = "ifany")))
saveRDS(pooled, .cache)
} # end cache

# =====================================================================
# SECTION 2 — DESCRIPTIVE TABLE (Supplementary Table 1)
# =====================================================================
analytic <- pooled |> filter(!is.na(edent), !is.na(byear), !is.na(male))
cat("\n== Analytic sample (complete edent/byear/male): N =", nrow(analytic), "==\n")

tab1 <- analytic |>
  mutate(age   = as.numeric(zap_labels(age)),
         byear = as.numeric(zap_labels(byear)),
         male  = as.numeric(zap_labels(male)),
         edent = as.numeric(zap_labels(edent)),
         coh = case_when(
           byear < famine_start - 5 ~ "Pre-famine",
           byear >= famine_start - 5 & byear <= famine_end ~ "Famine-exposed",
           byear > famine_end ~ "Post-famine"),
         coh = factor(coh, c("Pre-famine","Famine-exposed","Post-famine"))) |>
  group_by(study_f, coh) |>
  summarise(n = n(),
            age_m     = round(mean(age, na.rm = TRUE), 1),
            age_sd    = round(sd(age,   na.rm = TRUE), 1),
            byear_m   = round(mean(byear, na.rm = TRUE), 0),
            male_pct  = round(100*mean(male,  na.rm = TRUE), 1),
            edent_pct = round(100*mean(edent, na.rm = TRUE), 1),
            .groups = "drop")
names(tab1) <- c("study_f","cohort3_f","n","age_m","age_sd","byear_m","male_pct","edent_pct")
write.csv(tab1, paste0(TAB, "Table1_descriptives.csv"), row.names = FALSE)
cat("Descriptives grand total =", sum(tab1$n),
    "| range", min(tapply(tab1$n, tab1$study_f, sum)),
    "to", max(tapply(tab1$n, tab1$study_f, sum)), "\n")

# =====================================================================
# SECTION 3 — PRIMARY ANALYSIS: REGRESSION DISCONTINUITY IN TIME (RDiT)
# Hausman & Rapson (2018); rdrobust (Calonico et al. 2014);
# RDHonest (Armstrong & Kolesár 2018); Donut RDiT.
# =====================================================================
rdd <- analytic |> mutate(run = byear - (famine_start - 5))
cov_mat <- model.matrix(~ male + factor(father_ses) + factor(study), data = rdd)[, -1]

# -- 3.1 Main pooled RDiT ----
rdd_main <- rdrobust(rdd$edent, rdd$run, c = 0, p = 1, kernel = "triangular",
                     bwselect = "mserd", covs = cov_mat, masspoints = "adjust")
cat("\n== RDiT MAIN (pooled, MSE-RD, mass-points adjusted) ==\n"); summary(rdd_main)

# -- 3.2 Robustness specifications ----
rd_specs <- list(`CER-RD` = list(bw="cerrd",p=1,h=NULL),
                 `MSE p=2` = list(bw="mserd",p=2,h=NULL),
                 `h = 3 yr` = list(bw=NULL,p=1,h=3),
                 `h = 7 yr` = list(bw=NULL,p=1,h=7))
rdd_rob <- map_dfr(names(rd_specs), function(nm) {
  sp <- rd_specs[[nm]]
  args <- list(y=rdd$edent, x=rdd$run, c=0, p=sp$p, kernel="triangular",
               covs=cov_mat, masspoints="adjust")
  if (!is.null(sp$bw)) args$bwselect <- sp$bw
  if (!is.null(sp$h))  args$h <- sp$h
  r <- do.call(rdrobust, args)
  tibble(spec=nm, coef=r$coef["Conventional","Coeff"],
         se=r$se["Conventional","Std. Err."],
         lb=r$ci["Conventional","CI Lower"],
         ub=r$ci["Conventional","CI Upper"],
         p=r$pv["Conventional","P>|z|"])
})
donut_rows <- map_dfr(c(1,2), function(d) {
  rdd_d <- filter(rdd, run < -d | run > d)
  cov_d <- model.matrix(~ male + factor(father_ses) + factor(study), data=rdd_d)[,-1]
  fit <- tryCatch(rdrobust(rdd_d$edent, rdd_d$run, c=0, p=1, kernel="triangular",
                           bwselect="mserd", covs=cov_d, masspoints="adjust"),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  tibble(spec = sprintf("Donut RDiT (drop |run| <= %d)", d),
         coef=fit$coef["Conventional","Coeff"],
         se=fit$se["Conventional","Std. Err."],
         lb=fit$ci["Conventional","CI Lower"],
         ub=fit$ci["Conventional","CI Upper"],
         p=fit$pv["Conventional","P>|z|"])
})
rdd_rob <- bind_rows(rdd_rob, donut_rows)

rdh_fit <- tryCatch(RDHonest(edent ~ run, data = rdd[, c("edent","run")],
                             cutoff = 0, kern = "triangular",
                             opt.criterion = "MSE", sclass = "H"),
                    error = function(e) NULL)
if (!is.null(rdh_fit)) {
  rdd_rob <- bind_rows(rdd_rob, tibble(
    spec = "RDHonest (honest CI, M-estimated)",
    coef = rdh_fit$coefficients$estimate,
    se   = rdh_fit$coefficients$std.error,
    lb   = rdh_fit$coefficients$conf.low,
    ub   = rdh_fit$coefficients$conf.high,
    p    = NA_real_))
}

rdd_tab <- bind_rows(
  tibble(spec = "Main rdrobust (MSE-RD)",
         coef = rdd_main$coef["Conventional","Coeff"],
         se   = rdd_main$se["Conventional","Std. Err."],
         lb   = rdd_main$ci["Conventional","CI Lower"],
         ub   = rdd_main$ci["Conventional","CI Upper"],
         p    = rdd_main$pv["Conventional","P>|z|"]),
  rdd_rob)
write.csv(rdd_tab, paste0(TAB, "rdit_pooled_specs.csv"), row.names = FALSE)
cat("\n== RDiT SPECIFICATIONS (pooled) ==\n"); print(rdd_tab, digits = 4)

# -- 3.3 Study-specific RDiT estimates ----
rdd_studies <- map_dfr(1:10, function(s) {
  d <- filter(rdd, study == s); if (nrow(d) < 100) return(NULL)
  cs <- model.matrix(~ male, data = d)[, -1, drop = FALSE]
  tryCatch({
    r <- rdrobust(d$edent, d$run, c=0, p=1, kernel="triangular",
                  bwselect="mserd", covs=cs, masspoints="adjust")
    tibble(study = s, label = study_labels[s],
           est = r$coef["Conventional","Coeff"],
           lb  = r$ci["Conventional","CI Lower"],
           ub  = r$ci["Conventional","CI Upper"],
           p   = r$pv["Conventional","P>|z|"])
  }, error = function(e) NULL)
})
write.csv(rdd_studies, paste0(TAB, "rdd_by_study.csv"), row.names = FALSE)
cat("\n== RDiT BY STUDY ==\n"); print(rdd_studies)

# -- 3.4 Figure 1A: forest plot of study-specific + pooled RDiT ----
forest_rdd <- bind_rows(
  tibble(label = "Pooled (rdrobust)",
         est = rdd_main$coef["Conventional","Coeff"],
         lb  = rdd_main$ci["Conventional","CI Lower"],
         ub  = rdd_main$ci["Conventional","CI Upper"],
         p   = rdd_main$pv["Conventional","P>|z|"],
         is_pool = TRUE),
  rdd_studies |> mutate(is_pool = FALSE) |> select(label, est, lb, ub, p, is_pool)
) |> mutate(label_f = fct_inorder(label))

p_rdd <- ggplot(forest_rdd,
       aes(x = 100*est, y = fct_rev(label_f), xmin = 100*lb, xmax = 100*ub)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(linewidth = I(ifelse(is_pool, 1.0, 0.55))),
                 height = 0.22, colour = "black") +
  geom_point(aes(shape = ifelse(is_pool, 18L, 19L),
                 size  = I(ifelse(is_pool, 4.3, 2.5))),
             colour = "black", fill = "white") +
  scale_shape_identity() +
  labs(x = "Risk difference (pp), 95% CI", y = NULL,
       title = "Sharp RDD at famine-exposure cut-off — pooled and study-specific") +
  theme_pub(10) + theme(legend.position = "none")
ggsave(paste0(FIG, "Fig2_forest_RDD.png"), p_rdd, width = 7.2, height = 5.0, dpi = 300)

# -- 3.5 Figure 1B: RDiT specifications ----
p_rdit <- ggplot(rdd_tab,
       aes(x = 100*coef, y = fct_inorder(spec), xmin = 100*lb, xmax = 100*ub)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(size = 0.5) +
  labs(x = "Risk difference (pp), 95% CI", y = NULL,
       title = "RDiT — pooled, by specification",
       subtitle = "rdrobust (bias-corrected) + Donut RDiT + RDHonest (bias-aware)") +
  theme_pub(11)
ggsave(paste0(FIG, "Fig2b_RDiT_specs.png"), p_rdit, width = 8.5, height = 4.5, dpi = 300)

cat("\nSaved Fig2_forest_RDD.png, Fig2b_RDiT_specs.png\n")

# =====================================================================
# SECTION 4 — SECONDARY: REGION-OF-BIRTH × FAMINE-INTENSITY DiD
# CHARLS (Great Famine 1959-61) + CLHLS (Henan famine 1942-43)
# Province of birth, Chen & Zhou (2007) design.
# Unadjusted TWFE = primary; adjusted + BJS = sensitivity.
# =====================================================================

# -- Province codes & high-mortality provinces (CHARLS) ----
gb <- c('11'='Beijing','12'='Tianjin','13'='Hebei','14'='Shanxi',
        '15'='Inner Mongolia','21'='Liaoning','22'='Jilin','23'='Heilongjiang',
        '31'='Shanghai','32'='Jiangsu','33'='Zhejiang','34'='Anhui','35'='Fujian',
        '36'='Jiangxi','37'='Shandong','41'='Henan','42'='Hubei','43'='Hunan',
        '44'='Guangdong','45'='Guangxi','46'='Hainan','50'='Chongqing',
        '51'='Sichuan','52'='Guizhou','53'='Yunnan','54'='Tibet',
        '61'='Shaanxi','62'='Gansu','63'='Qinghai','64'='Ningxia','65'='Xinjiang')
HIGH <- c('Anhui','Sichuan','Chongqing','Guizhou','Gansu','Qinghai',
          'Henan','Hunan','Guangxi')

# -- CHARLS: extract birth province ----
res <- read_dta(paste0(DAT, "CHARLS/2014_life_history/Residence.dta"))
res$psu1 <- as.numeric(res$t002_psu_1_)
bp <- tibble(ID = res$ID,
             birthprov = gb[as.character(ifelse(res$psu1 >= 1100,
                                               floor(res$psu1/100), NA))])

charls_did <- read_dta(paste0(DAT, "CHARLS/harmonized/H_CHARLS_D_Data.dta")) |>
  transmute(byear = rabyear, age = r1agey,
            male = if_else(ragender == 1, 1L, 0L), ID) |>
  left_join(read_dta(paste0(DAT, "CHARLS/w2/Health_Status_and_Functioning.dta")) |>
              transmute(ID, edent = case_when(da040 == 1 ~ 1L, da040 == 2 ~ 0L)), "ID") |>
  left_join(bp, "ID") |>
  transmute(study = "CHARLS", famine = "Great Famine 1959–61",
            edent, byear, age, male, region = birthprov,
            treat = as.integer(birthprov %in% HIGH),
            gstart = 1954L, famine_lo = 1959L, famine_hi = 1961L,
            pre = byear %in% 1947:1953,
            exposed = byear %in% 1954:1960,
            post = byear %in% 1961:1967)

# -- CLHLS: extract birth province ----
clhls_did <- read_dta(paste0(DAT,
  "CLHLS/clhls_2008_2018_longitudinal_dataset_released_version1.dta")) |>
  transmute(study = "CLHLS", famine = "Henan famine 1942–43",
            edent = as.integer(g21 == 0 & !is.na(g21) & g21 <= 32),
            byear = v_bthyr, age = trueage, male = as.integer(a1 == 1),
            region = as.character(as_factor(a41))) |>
  mutate(treat = as.integer(region == "henan"),
         gstart = 1937L, famine_lo = 1942L, famine_hi = 1943L,
         pre = byear %in% 1931:1936,
         exposed = byear %in% 1937:1942,
         post = byear %in% 1943:1948)

studies_did <- list(CHARLS = charls_did, CLHLS = clhls_did)

# -- TWFE estimator (unadjusted + adjusted) ----
twfe <- function(d, adj = TRUE) {
  d <- d[complete.cases(d[, c("edent","age","male","region","treat")]), ]
  d$exposed <- as.integer(d$exposed)
  f <- if (adj) lm(edent ~ treat*exposed + age + I(age^2) + I(age^3) + male, data = d)
       else     lm(edent ~ treat*exposed, data = d)
  coeftest(f, vcovCL(f, cluster = d$region))["treat:exposed", c(1,2,4)]
}
bjs <- function(s) {
  sb <- s |> mutate(Gi = ifelse(treat == 1, gstart[1], 0))
  b <- did_imputation(data = sb, yname = "edent", gname = "Gi", tname = "byear",
                      first_stage = ~ male | region + byear, cluster_var = "region")
  c(b$estimate[1], b$std.error[1], 2*pnorm(-abs(b$estimate[1]/b$std.error[1])))
}

est <- list(); es_all <- list()
for (nm in names(studies_did)) {
  d <- studies_did[[nm]]
  cat(sprintf("\n== Region-of-birth DiD: %s ==\n", nm))
  s <- d |> filter(pre | exposed, !is.na(edent), !is.na(age), !is.na(male), !is.na(region))
  rU <- twfe(s, adj = FALSE); rT <- twfe(s, adj = TRUE)
  rP <- twfe(d |> filter(exposed | post, !is.na(edent), !is.na(age), !is.na(male), !is.na(region)),
             adj = TRUE)
  rB <- bjs(s)
  ci <- function(r) sprintf("%+5.2f to %+5.2f", 100*(r[1]-1.96*r[2]), 100*(r[1]+1.96*r[2]))
  cat(sprintf("  TWFE unadjusted (primary) %+6.2f pp (95%% CI %s, p=%.3f)\n", 100*rU[1], ci(rU), rU[3]))
  cat(sprintf("  TWFE adjusted (sensitivity) %+6.2f pp (95%% CI %s, p=%.3f)\n", 100*rT[1], ci(rT), rT[3]))
  cat(sprintf("  TWFE post control (sensitivity) %+6.2f pp (SE %.2f, p=%.3f)\n", 100*rP[1], 100*rP[2], rP[3]))
  cat(sprintf("  BJS imputation %+6.2f pp (SE %.2f, p=%.3f)\n", 100*rB[1], 100*rB[2], rB[3]))
  est[[nm]] <- rbind(
    data.frame(study=nm, estimator="TWFE unadjusted (primary)", est=100*rU[1], se=100*rU[2], p=rU[3]),
    data.frame(study=nm, estimator="TWFE adjusted (sensitivity)", est=100*rT[1], se=100*rT[2], p=rT[3]),
    data.frame(study=nm, estimator="BJS imputation (modern)", est=100*rB[1], se=100*rB[2], p=rB[3]),
    data.frame(study=nm, estimator="TWFE post control (sensitivity)", est=100*rP[1], se=100*rP[2], p=rP[3]))

  # Event study (cohort-by-region interaction)
  gs <- d$gstart[1]; refyr <- gs - 1L
  esd <- d |> filter(byear >= gs-7 & byear <= gs+13,
                     !is.na(edent), !is.na(male), !is.na(region))
  m <- feols(edent ~ i(byear, treat, ref = refyr) + male | region + byear,
             data = esd, cluster = ~region)
  ct <- as.data.frame(coeftable(m)); ct <- ct[grepl("byear::", rownames(ct)), ]
  ct$byear <- as.integer(sub(":treat", "", sub("byear::", "", rownames(ct))))
  ct$rel <- ct$byear - gs; ct$study <- nm; ct$famine <- d$famine[1]
  ct$famine_lo <- d$famine_lo[1]; ct$famine_hi <- d$famine_hi[1]
  names(ct)[1:2] <- c("coef","se")
  es_all[[nm]] <- ct
}
estdf <- do.call(rbind, est)
write.csv(estdf, paste0(TAB, "region_did_final.csv"), row.names = FALSE)
esdf <- do.call(rbind, es_all)
write.csv(esdf[, c("study","famine","byear","rel","coef","se")],
          paste0(TAB, "region_eventstudy.csv"), row.names = FALSE)
cat("\n== Region-DiD estimates ==\n"); print(estdf)

# -- Figure 2: event-study panels with famine band ----
esdf$lab <- paste0(esdf$study, " — ", esdf$famine)
band <- do.call(rbind, lapply(split(esdf, esdf$lab), function(z)
  data.frame(lab = z$lab[1], xmin = z$famine_lo[1]-0.5, xmax = z$famine_hi[1]+0.5)))
cutl <- do.call(rbind, lapply(split(esdf, esdf$lab), function(z)
  data.frame(lab = z$lab[1], gs = (z$byear[1]-z$rel[1])-0.5)))

p_did <- ggplot(esdf, aes(byear, 100*coef)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.12) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(data = cutl, aes(xintercept = gs), linetype = 3, colour = "firebrick") +
  geom_pointrange(aes(ymin = 100*(coef-1.96*se), ymax = 100*(coef+1.96*se)), size = .3) +
  facet_wrap(~lab, scales = "free_x", ncol = 1) +
  labs(x = "Birth year", y = "Treated−control difference in edentulism (pp)") +
  theme_pub(10)
ggsave(paste0(FIG, "Fig_region_did.png"), p_did, width = 7.0, height = 5.0, dpi = 300)
cat("\nSaved Fig_region_did.png\n")

# =====================================================================
# SECTION 5 — SECONDARY: INTERRUPTED TIME-SERIES (ITA)
# Segmented regression with level shift + slope change at famine cohort.
# =====================================================================
its <- pooled |>
  filter(!is.na(edent), !is.na(byear)) |>
  group_by(study, byear, famine_start, famine_end) |>
  summarise(edent = mean(edent, na.rm = TRUE), n = sum(!is.na(edent)), .groups = "drop") |>
  mutate(time = byear - famine_start,
         post = as.integer(byear >= famine_start - 5),
         time_post = time * post)
ITS_mod <- lm(edent ~ time + post + time_post + factor(study), data = its, weights = n)
its_out <- as.data.frame(coeftest(ITS_mod, vcov = vcovHC(ITS_mod, type = "HC3"))[1:4, ])
write.csv(its_out, paste0(TAB, "ita_pooled.csv"), row.names = TRUE)
cat("\n== ITA (level shift + slope change at famine-exposure cut-off) ==\n")
print(its_out, digits = 4)

# -- Figure 3: ITA by study ----
its_p <- its |>
  mutate(event_time = byear - famine_start) |>
  filter(event_time >= -20, event_time <= 20)
p_ita <- ggplot(its_p, aes(event_time, 100*edent)) +
  geom_vline(xintercept = -5, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = n), alpha = 0.35) +
  geom_smooth(data = ~ filter(.x, event_time < -5),
              method = "lm", se = FALSE, colour = "#2C7FB8") +
  geom_smooth(data = ~ filter(.x, event_time >= -5),
              method = "lm", se = FALSE, colour = "#E69F00") +
  facet_wrap(~ study, scales = "free_y",
             labeller = labeller(study = setNames(study_labels, as.character(1:10)))) +
  scale_size_continuous(range = c(0.3, 2.5), guide = "none") +
  labs(x = "Event time (birth year − famine onset)",
       y = "Edentulism prevalence (%)",
       title = "Interrupted Time-Series — birth-year trends by study") +
  theme_pub(9)
ggsave(paste0(FIG, "Fig4_ITA.png"), p_ita, width = 9, height = 6.5, dpi = 300)
cat("\nSaved Fig4_ITA.png\n")

# =====================================================================
# SECTION 6 — NEGATIVE CONTROLS
# (a) Placebo cohort (-10 yrs), (b) Placebo window (age 11-15),
# (c) Negative-control outcome (sex / not affected by famine).
# =====================================================================
nc1 <- pooled |>
  mutate(p_start = famine_start-10, p_end = famine_end-10,
         p_treat = case_when(byear >= p_start-5 & byear <= p_end-1 ~ 1L,
                             byear >= p_end+1 & byear <= p_end+6 ~ 0L)) |>
  filter(!is.na(p_treat), !is.na(edent))
nc1_row <- tidy(feols(edent ~ p_treat + age_c + age_c2 + age_c3 + male + i(father_ses) | study_f,
                      data = nc1, cluster = ~study_f)) |>
  filter(term == "p_treat") |> mutate(spec = "Placebo cohort (-10 yrs)")

nc2 <- pooled |>
  mutate(p2_treat = case_when(byear >= famine_start-15 & byear <= famine_end-11 ~ 1L,
                              byear >= famine_start-20 & byear <= famine_start-16 ~ 0L)) |>
  filter(!is.na(p2_treat), !is.na(edent))
nc2_row <- tidy(feols(edent ~ p2_treat + age_c + age_c2 + age_c3 + male + i(father_ses) | study_f,
                      data = nc2, cluster = ~study_f)) |>
  filter(term == "p2_treat") |> mutate(spec = "Placebo window (age 11-15)")

nc3_row <- tidy(feols(male ~ treatment + age_c + age_c2 + age_c3 + i(father_ses) | study_f,
                      data = filter(pooled, !is.na(treatment), !is.na(male)),
                      cluster = ~study_f)) |>
  filter(term == "treatment") |> mutate(spec = "Neg.-ctrl outcome (male)")

nc_tab <- bind_rows(nc1_row, nc2_row, nc3_row) |>
  mutate(est_pp = 100*estimate,
         lb_pp  = 100*(estimate - 1.96*std.error),
         ub_pp  = 100*(estimate + 1.96*std.error)) |>
  select(spec, estimate, std.error, p.value, est_pp, lb_pp, ub_pp)
write.csv(nc_tab, paste0(TAB, "negative_controls.csv"), row.names = FALSE)
cat("\n== NEGATIVE CONTROLS ==\n"); print(nc_tab, digits = 3)

# =====================================================================
# SECTION 7 — SUPPLEMENT SENSITIVITY ANALYSES
# S2: Reduced-form RR attenuation; S3: Leave-one-out; S4: Duration.
# =====================================================================

# -- S2: Reduced-form models (risk differences via LPM) ----
dat2 <- pooled |> filter(cohort3 %in% 1:2)
M_orig  <- lm(edent ~ factor(treatment) + factor(study), data = dat2)
M_adj   <- lm(edent ~ factor(treatment) + male + factor(father_ses) + factor(study), data = dat2)
M_ageaj <- lm(edent ~ factor(treatment) + age_c + age_c2 + age_c3 +
               male + factor(father_ses) + factor(study), data = dat2)
rf_supp <- map_dfr(list(Unadjusted = M_orig, `+ Sex + SES` = M_adj,
                        `+ Age cubic` = M_ageaj),
                   function(m) {
                     r  <- coef(summary(m))["factor(treatment)1", ]
                     ci <- confint(m)["factor(treatment)1", ]
                     tibble(coef=r[1], se=r[2], lb=ci[1], ub=ci[2], p=r[4], n=nobs(m))
                   }, .id = "Model")
write.csv(rf_supp, paste0(TAB, "S1_reduced_form_models.csv"), row.names = FALSE)
cat("\n== SUPPLEMENT S2 — Reduced-form models ==\n"); print(rf_supp, digits = 4)

# -- S3: Leave-one-out ----
loo <- map_dfr(1:10, function(s) {
  d <- filter(dat2, study != s)
  m <- lm(edent ~ factor(treatment) + age_c + age_c2 + age_c3 + male +
           factor(father_ses) + factor(study), data = d)
  r  <- coef(summary(m))["factor(treatment)1", ]
  ci <- confint(m)["factor(treatment)1", ]
  tibble(Excluded = study_labels[s], coef=r[1], se=r[2],
         lb=ci[1], ub=ci[2], p=r[4], n=nobs(m))
})
write.csv(loo, paste0(TAB, "S2_leave_one_out.csv"), row.names = FALSE)
cat("\n== SUPPLEMENT S3 — Leave-one-out ==\n"); print(loo, digits = 4)

# -- S4: Famine-duration stratification ----
short_s <- pooled |> distinct(study, short_famine) |> filter(short_famine == 1) |> pull(study)
long_s  <- pooled |> distinct(study, short_famine) |> filter(short_famine == 0) |> pull(study)
dur <- map_dfr(list(Short = short_s, Long = long_s), function(ss) {
  d <- filter(dat2, study %in% ss)
  m <- lm(edent ~ factor(treatment) + age_c + age_c2 + age_c3 + male +
           factor(father_ses) + factor(study), data = d)
  r  <- coef(summary(m))["factor(treatment)1", ]
  ci <- confint(m)["factor(treatment)1", ]
  tibble(coef=r[1], se=r[2], lb=ci[1], ub=ci[2], p=r[4], n=nobs(m))
}, .id = "Duration")
write.csv(dur, paste0(TAB, "S3_duration_strat.csv"), row.names = FALSE)
cat("\n== SUPPLEMENT S4 — Duration stratification ==\n"); print(dur, digits = 4)

# =====================================================================
# SESSION INFO
# =====================================================================
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\n=================== ANALYSIS COMPLETE ===================\n")
cat("Tables   -> tables/\n")
cat("Figures  -> figures/\n")
cat("Session  -> sessionInfo.txt\n")
