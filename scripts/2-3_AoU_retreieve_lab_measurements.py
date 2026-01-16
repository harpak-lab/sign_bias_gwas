## Get WBC lab measurements (monocyte, basophil, & neutrophil percentages)

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dataset = %env WORKSPACE_CDR
bucket = os.getenv("WORKSPACE_BUCKET")

#make df
lab_df = pd.read_gbq(f'''SELECT DISTINCT
  m.person_id,
  LOWER(m.measurement_source_value) AS measurement_source_value,
  m.value_as_number,
  COALESCE(u.concept_name, CAST(m.unit_source_value AS STRING)) AS unit_concept_name
FROM `{dataset}.measurement` AS m
LEFT JOIN `{dataset}.concept` AS u
  ON u.concept_id = m.unit_concept_id
WHERE LOWER(m.measurement_source_value) IN ('706-2','770-8','5905-5')''')

loinc_to_name = {
    '706-2': 'basophil_percentage',
    '770-8': 'neutrophil_percentage',
    '5905-5': 'monocyte_percentage',
}
lab_df['lab'] = lab_df['measurement_source_value'].map(loinc_to_name)

lab_df['unit_concept_name'].astype(str).str.lower().value_counts()

#remove invalid concept names
keep = ['percent', 'percentage unit', 'percent of white blood cells']
lab_df_n = lab_df[lab_df['unit_concept_name'].isin(keep)]
lab_df_b = lab_df_n.copy()

# Default bounds 
bounds = {
    'basophil_percentage':  (0, 100),  
    'neutrophil_percentage':(0, 100),  
    'monocyte_percentage':  (0, 100),   
}

# Helper to change bounds later:
def set_bounds(new_bounds: dict):
    bounds.update(new_bounds)

# Attach bounds row-wise
lab_df_b[['lb','ub']] = pd.DataFrame(
    lab_df_b['lab'].map(lambda k: bounds.get(k, (0,100))).tolist(),
    index=lab_df_b.index
)

lab_df_b['is_outlier'] = (lab_df_b['value_as_number'] < lab_df_b['lb']) | (lab_df_b['value_as_number'] > lab_df_b['ub'])
lab_df_b['is_missing'] = lab_df_b['value_as_number'].isna()

clean = lab_df_b.loc[~lab_df_b['is_missing'] & ~lab_df_b['is_outlier']].copy()

summary = (clean
           .groupby(['person_id','lab'])['value_as_number']
           .agg(min='min', median='median', max='max', mean='mean', count='size')
           .reset_index())

summary.to_csv("lab_measures_wbc.tsv", sep = "\t")

#copy to bucket
!gsutil -m -u $GOOGLE_PROJECT cp lab_measures_wbc.tsv '{bucket}/aou_gwas/pheno/'