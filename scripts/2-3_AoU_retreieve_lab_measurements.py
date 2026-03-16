## Get lab measurements
###########

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
WHERE LOWER(m.measurement_source_value) IN ('706-2','770-8','5905-5','6690-2','789-8', '785-6')''')

#rename from loinc
loinc_to_name = {
    '706-2': 'basophil_percentage',
    '770-8': 'neutrophil_percentage',
    '5905-5': 'monocyte_percentage',
    '6690-2': 'white_blood_cell_count',
    '789-8': 'red_blood_cell_count',
    '785-6': "mean_corpuscular_hemoglobin"
}
lab_df['lab'] = lab_df['measurement_source_value'].map(loinc_to_name)

# get raw stats for each trait
raw_stats = (lab_df
             .groupby('lab', dropna=False)['value_as_number']
             .agg(n='size', n_missing=lambda x: x.isna().sum(),
                  min='min', median='median', max='max')
             .sort_values('max', ascending=False))

#percentage unit names
percent_labs = {'basophil_percentage','neutrophil_percentage','monocyte_percentage'}
percent_units_keep = {'percent', 'percentage unit', 'percent of white blood cells'}  

df_pct = lab_df[lab_df['lab'].isin(percent_labs)].copy()
df_pct = df_pct[df_pct['unit_concept_name'].isin(percent_units_keep)].copy()
df_pct = df_pct[df_pct['value_as_number'].notna()].copy()

#bound percentages 0-100
percent_labs = ['basophil_percentage','neutrophil_percentage','monocyte_percentage']
percent_units_keep = ['percent', 'percentage unit', 'percent of white blood cells']

df_pct = lab_df[lab_df['lab'].isin(percent_labs)].copy()
df_pct = df_pct[df_pct['unit_concept_name'].isin(percent_units_keep)].copy()

df_pct['value_std'] = df_pct['value_as_number']
df_pct['unit_std'] = 'percent'

bounds = {lab: (0,100) for lab in percent_labs}

df_pct[['lb','ub']] = pd.DataFrame(
    df_pct['lab'].map(lambda k: bounds.get(k, (0,100))).tolist(),
    index=df_pct.index
)

df_pct['is_missing'] = df_pct['value_std'].isna()
df_pct['is_outlier'] = (df_pct['value_std'] < df_pct['lb']) | (df_pct['value_std'] > df_pct['ub'])

clean_pct = df_pct.loc[~df_pct['is_missing'] & ~df_pct['is_outlier']].copy()

#summary of percentage traits
summary_pct = (clean_pct
               .groupby(['person_id','lab'])['value_as_number']
               .agg(min='min', median='median', max='max', mean='mean', count='size')
               .reset_index())
summary_pct.head()

#count unit names for blood cell counts
count_labs = {'white_blood_cell_count','red_blood_cell_count'}

# make unit names consistent
df_cnt = lab_df[lab_df['lab'].isin(count_labs)].copy()
df_cnt['unit_concept_name'] = df_cnt['unit_concept_name'].astype(str).str.lower().str.strip()
df_cnt = df_cnt[df_cnt['value_as_number'].notna()].copy()

# unit harmonization
count_unit_maps = {
    "white_blood_cell_count": {
        "unit_std": "thousand per microliter",
        "unit_factors": {
            "thousand per microliter": 1.0,
            "thousand per cubic millimeter": 1.0,   
            "x10(3)/mcl": 1.0,
            "billion per liter": 1.0,               

            "cells per microliter": 1e-3,          
            "cells/ul": 1e-3,
            "per microliter": 1e-3,
            "per cubic millimeter": 1e-3,
            "/mm3": 1e-3,
        }
    },
    "red_blood_cell_count": {
        "unit_std": "million per microliter",
        "unit_factors": {
            "million per microliter": 1.0,
            "million per cubic millimeter": 1.0,

            "cells per microliter": 1e-6,           
            "cells/ul": 1e-6,
            "per microliter": 1e-6,
            "per cubic millimeter": 1e-6,
            "/mm3": 1e-6,

            "thousand per microliter": 1e-3,      
            "thousand per cubic millimeter": 1e-3,
            "billion per liter": 1e-3,             
            "million per liter": 1e-6,            
        }
    }
}

# mapping table
map_rows = []
for lab, spec in count_unit_maps.items():
    for u, f in spec["unit_factors"].items():
        map_rows.append({"lab": lab, "unit_concept_name": u, "factor": f, "unit_std": spec["unit_std"]})
count_map = pd.DataFrame(map_rows)

# keep 
df_cnt = df_cnt.merge(count_map, on=["lab","unit_concept_name"], how="left")
df_cnt = df_cnt[df_cnt["factor"].notna()].copy()

# standardize
df_cnt["value_std"] = df_cnt["value_as_number"] * df_cnt["factor"]

#bound blood cell counts
abs_bounds = {
    "white_blood_cell_count": (0.1, 500.0),  
    "red_blood_cell_count":   (0.5, 12.0),   
}

df_cnt[["lb","ub"]] = pd.DataFrame(
    df_cnt["lab"].map(lambda k: abs_bounds.get(k, (np.nan, np.nan))).tolist(),
    index=df_cnt.index
)

df_cnt["is_missing"] = df_cnt["value_std"].isna()
df_cnt["is_outlier"] = (df_cnt["value_std"] < df_cnt["lb"]) | (df_cnt["value_std"] > df_cnt["ub"]) | (df_cnt["value_std"] <= 0)

clean_cnt = df_cnt.loc[~df_cnt["is_missing"] & ~df_cnt["is_outlier"]].copy()

#summarize counts
summary_cnt = (clean_cnt
               .groupby(["person_id","lab"])["value_std"]
               .agg(min="min", median="median", max="max", mean="mean", count="size")
               .reset_index())

summary_cnt.head()

#keep mch units
mch_units_keep = {'picogram', 'picogram per cell', 'pg/cell'}

df_mch = lab_df[lab_df['lab'] == 'mean_corpuscular_hemoglobin'].copy()
df_mch = df_mch[df_mch['unit_concept_name'].isin(mch_units_keep)].copy()

df_mch['value_std'] = df_mch['value_as_number']
df_mch['unit_std']  = 'pg'

#bound units
bounds = {'mean_corpuscular_hemoglobin': (5, 100)}

df_mch[['lb','ub']] = pd.DataFrame(
    df_mch['lab'].map(lambda k: bounds.get(k, (5,80))).tolist(),
    index=df_mch.index
)

df_mch['is_missing'] = df_mch['value_std'].isna()
df_mch['is_outlier'] = (df_mch['value_std'] < df_mch['lb']) | (df_mch['value_std'] > df_mch['ub'])

clean_mch = df_mch.loc[~df_mch['is_missing'] & ~df_mch['is_outlier']].copy()

mch_median = (clean_mch
              .groupby('person_id')['value_std']
              .median()
              .reset_index(name='mean_corpuscular_hemoglobin_median'))


# percentage traits: wide format
pct_wide = (
    summary_pct[['person_id', 'lab', 'median']]
    .pivot(index='person_id', columns='lab', values='median')
    .reset_index()
)

# rename cols
pct_wide = pct_wide.rename(columns={
    'basophil_percentage': 'basophil_percentage_median',
    'monocyte_percentage': 'monocyte_percentage_median',
    'neutrophil_percentage': 'neutrophil_percentage_median'
})

# count traits: wide format
cnt_wide = (
    summary_cnt[['person_id', 'lab', 'median']]
    .pivot(index='person_id', columns='lab', values='median')
    .reset_index()
)

# rename cols
cnt_wide = cnt_wide.rename(columns={
    'white_blood_cell_count': 'white_blood_cell_count_median',
    'red_blood_cell_count': 'red_blood_cell_count_median'
})

# merge 
final_labs = (
    pct_wide
    .merge(cnt_wide, on='person_id', how='outer')
    .merge(mch_median, on='person_id', how='outer')
)

# order
final_labs = final_labs[[
    'person_id',
    'basophil_percentage_median',
    'mean_corpuscular_hemoglobin_median',
    'monocyte_percentage_median',
    'neutrophil_percentage_median',
    'red_blood_cell_count_median',
    'white_blood_cell_count_median'
]]

#write out
final_labs.to_csv(
    "lab_measures_all_v2.tsv",
    sep="\t",
    index=False,
    na_rep="NA"
)

#copy to bucket
!gsutil -m -u $GOOGLE_PROJECT cp lab_measures_all_v2.tsv '{bucket}/aou_gwas/pheno/'