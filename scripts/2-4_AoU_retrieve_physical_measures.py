## Get AoU physical measurements (height, weight, and BMI)

import os
import numpy as np
import pandas as pd
dataset = %env WORKSPACE_CDR

calculated_pm_df = pd.read_gbq(f'''SELECT DISTINCT
  m.person_id,
  LOWER(m.measurement_source_value) AS measurement_source_value,
  m.value_as_number,
  COALESCE(u.concept_name, CAST(m.unit_source_value AS STRING)) AS unit_concept_name
FROM `{dataset}.measurement` AS m
LEFT JOIN `{dataset}.concept` AS u
  ON u.concept_id = m.unit_concept_id
WHERE LOWER(m.measurement_source_value) IN ('weight','height')''')

# remove rows with no matching concept in unit_concept_name
calculated_pm_df = calculated_pm_df[calculated_pm_df.unit_concept_name!="No matching concept"]

pm_df = (
    calculated_pm_df
      .drop(columns='unit_concept_name')
      .pivot_table(index='person_id',
                   columns='measurement_source_value',
                   values='value_as_number')
      .reset_index()
      .rename_axis(None, axis=1)        
      .rename(columns={'height': 'height_cm', 'weight': 'weight_kg'})
)

# Get BMI
# BMI = weight_kg / (height_m^2); height_m = height_cm / 100
pm_df["bmi"] = np.where(
    pm_df["height_cm"] > 0,
    pm_df["weight_kg"] / (pm_df["height_cm"] / 100.0) ** 2,
    np.nan
)

pm_df[["height_cm", "weight_kg", "bmi"]].describe()

pm_df.to_csv("physical_measures_hwb.tsv", sep = "\t")

## Copy to bucket
bucket = os.getenv("WORKSPACE_BUCKET")
!gsutil -m -u $GOOGLE_PROJECT cp physical_measures_hwb.tsv '{bucket}/aou_gwas/pheno/'
