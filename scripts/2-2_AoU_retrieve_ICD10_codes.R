#### Get ICD10 codes from PhecodeX tables

library(reticulate)

use_python("/opt/conda/bin/python", required = TRUE)
reticulate::py_discover_config()
check <- py_capture_output({
  py_run_string("import sys, ssl; print(sys.executable); print(ssl.OPENSSL_VERSION)")
})
cat(check)

# load libraries
library(tidyverse)
if (!requireNamespace("PheWAS", quietly = TRUE)) devtools::install_github("PheWAS/PheWAS")
library(PheWAS)
library(bigrquery)
library(glue)
library(readr)

#Set up workspace

BILLING_PROJECT_ID <- Sys.getenv("GOOGLE_PROJECT")
CDR                <- Sys.getenv("WORKSPACE_CDR")
WORKSPACE_BUCKET   <- Sys.getenv("WORKSPACE_BUCKET")
git_base <- "https://raw.githubusercontent.com/PheWAS/PhecodeXVocabulary/main/PhecodeX%20(version%201.0)"
map_url  <- glue("{git_base}/phecodeX_unrolled_ICD_CM.csv")
map_df   <- readr::read_csv(map_url, show_col_types = FALSE)
head(map_df)

readr::write_csv(map_df, "updated_phecodex_map.csv")
system2("gsutil", args = c("cp", "updated_phecodex_map.csv",
                           paste0(WORKSPACE_BUCKET, "/aou_gwas/phecodeX_files/")),
        stdout = TRUE, stderr = TRUE)

bigquery <- reticulate::import("google.cloud.bigquery")
client   <- bigquery$Client()
phecode <- bigquery$ExternalConfig("CSV")
phecode$source_uris <- paste0(WORKSPACE_BUCKET, "/aou_gwas/phecodeX_files/updated_phecodex_map.csv")
phecode$schema <- reticulate::r_to_py(list(
  bigquery$SchemaField("phecode", "STRING"),
  bigquery$SchemaField("ICD",  "STRING"),
  bigquery$SchemaField("vocabulary_id",       "STRING")
))
reticulate::py_set_attr(phecode$options, "skip_leading_rows", "1")

#Query

job_config <- bigquery$QueryJobConfig()
job_config$default_dataset   <- CDR
job_config$table_definitions <- reticulate::dict("phecodex" = phecode)
peek <- client$query("SELECT vocabulary_id, ICD, phecode FROM phecodex LIMIT 5",
                     job_config = job_config)$to_dataframe()


icds <- "
WITH all_codes AS (
  SELECT DISTINCT
    co.person_id,
    c.vocabulary_id,
    c.concept_code,
    co.condition_start_date AS date
  FROM `condition_occurrence` co
  JOIN `concept` c
    ON c.concept_id = co.condition_source_concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')

  UNION DISTINCT

  SELECT DISTINCT
    o.person_id,
    c.vocabulary_id,
    c.concept_code,
    o.observation_date AS date
  FROM `observation` o
  JOIN `concept` c
    ON c.concept_id = o.observation_source_concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')

  UNION DISTINCT

  SELECT DISTINCT
    p.person_id,
    c.vocabulary_id,
    c.concept_code,
    p.procedure_date AS date
  FROM `procedure_occurrence` p
  JOIN `concept` c
    ON c.concept_id = p.procedure_source_concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')

  UNION DISTINCT

  SELECT DISTINCT
    m.person_id,
    c.vocabulary_id,
    c.concept_code,
    m.measurement_date AS date
  FROM `measurement` m
  JOIN `concept` c
    ON c.concept_id = m.measurement_source_concept_id
  WHERE c.vocabulary_id IN ('ICD9CM','ICD10CM')
)

SELECT
  a.person_id,
  p.phecode,
  COUNT(DISTINCT a.date) AS code_count
FROM all_codes a
JOIN (
  SELECT
    -- your external CSV has columns: vocabulary_id, ICD, phecode
    REGEXP_REPLACE(UPPER(vocabulary_id), r'[^A-Z0-9]', '') AS vocabulary_id,
    REPLACE(UPPER(ICD), '.', '')                           AS concept_code,
    phecode
  FROM phecodex
) AS p
  ON REGEXP_REPLACE(UPPER(a.vocabulary_id), r'[^A-Z0-9]', '') = p.vocabulary_id
 AND REPLACE(UPPER(a.concept_code), '.', '')                  = p.concept_code
GROUP BY a.person_id, p.phecode
"
job_config <- job_config  # your existing job_config with the external CSV
query_job <- client$query(icds, job_config = job_config)
data <- query_job$to_dataframe()
nrow(data)

job_config2 <- bigquery$QueryJobConfig(); job_config2$default_dataset <- CDR
ehr_qry <- client$query("SELECT DISTINCT person_id FROM cb_search_person WHERE has_ehr_data = 1",
                        job_config = job_config2)
ehr_inds <- ehr_qry$to_dataframe()
ehr_person_ids <- unique(ehr_inds$person_id)

has_count_col <- "count.column" %in% names(formals(PheWAS::createPhenotypes))
phe_input <- dplyr::as_tibble(data) %>%
  dplyr::transmute(
    person_id,
    vocabulary_id = "phecode",  
    phecode,
    code_count    = as.integer(code_count)
  )
head(phe_input)

# Call createPhenotypes 
phe_table <- createPhenotypes(
  phe_input,
  translate             = FALSE,
  min.code.count        = 2,
  add.phecode.exclusions= FALSE,
  full.population.ids   = ehr_person_ids
)

#  Convert to explicit booleans (TRUE/FALSE) 
phe_bool <- phe_table %>%
  rename(person_id = 1) %>%
  mutate(across(-person_id, ~ {
    x <- .
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) x <- tolower(x)
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x == 1)
    (x %in% c("1","case","true"))
  }))

# Sex column
sex_try <- try({
  sex_q <- "
    SELECT DISTINCT person_id,
      CASE
        WHEN UPPER(sex_at_birth) LIKE 'F%' THEN 'F'
        WHEN UPPER(sex_at_birth) LIKE 'M%' THEN 'M'
        ELSE NULL
      END AS sex
    FROM cb_search_person
  "
  as_tibble(client$query(sex_q, job_config = job_config2)$to_dataframe())
}, silent = TRUE)

if (!inherits(sex_try, "try-error") && all(c("person_id","sex") %in% names(sex_try))) {
  phe_bool <- phe_bool %>% left_join(sex_try, by = "person_id")
}

# Order columns: person_id, phecodes sorted, sex
phe_cols <- setdiff(names(phe_bool), c("person_id","sex"))
phe_bool <- phe_bool %>% select(person_id, all_of(sort(phe_cols)), any_of("sex"))

# Write output (local). 
readr::write_csv(phe_bool, "mcc2_phecodex_table_v8.csv")
message("Wrote mcc2_phecodex_table_v8.csv with ", nrow(phe_bool), " participants and ",
        length(phe_cols), " phecodeX columns.")