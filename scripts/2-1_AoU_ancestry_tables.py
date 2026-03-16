## Retrieve ancestry information and PCs for AoU participants

import pandas as pd
dataset = %env WORKSPACE_CDR
demographics_df = pd.read_gbq(f'''SELECT DISTINCT person_id, age_at_cdr, sex_at_birth, race 
                                FROM `{dataset}.cb_search_person`''')
age_bins = [0,29,39,49,59,69,79,89,200]
age_labels = ["18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"]
demographics_df['age_group'] = pd.cut(demographics_df['age_at_cdr']
                                      , bins=age_bins, labels=age_labels, include_lowest=True)

demographics_df.to_csv('demographics_table.csv', index = False)
genomics_ancestry_bucket = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry'
!gsutil -u $GOOGLE_PROJECT cp '{genomics_ancestry_bucket}/ancestry_preds.tsv' './ancestry.tsv'
ancestry_table = pd.read_csv('ancestry.tsv',sep='\t')
ancestry_table = ancestry_table[['research_id','ancestry_pred']]
ancestry_table.columns = ['person_id', 'ancestry']
ancestry_table
ancestry_table.to_csv('ancestry_table.csv', index = False)
ancestry_table = pd.read_csv('ancestry.tsv',sep='\t')
ancestry_table
import ast, numpy as np

def to_seq(x):
    if isinstance(x, (list, tuple, np.ndarray)): return list(x)
    if pd.isna(x): return []
    try: return list(ast.literal_eval(x))  # parse stringified list
    except Exception: return []

seqs = ancestry_table['pca_features'].map(to_seq)
max_len = seqs.map(len).max()
pca_cols = [f'PC{i+1}' for i in range(max_len)]

# pad shorter rows with NaN so they’re all equal length
pca_df = pd.DataFrame(
    seqs.apply(lambda xs: xs + [np.nan]*(max_len - len(xs))).tolist(),
    columns=pca_cols,
    index=ancestry_table.index
)

ancestry_with_pc = ancestry_table.drop(columns=['pca_features']).join(pca_df)
ancestry_with_pc = ancestry_with_pc.rename(columns={ancestry_with_pc.columns[0]: 'person_id'})
ancestry_with_pc.to_csv('participant_PCs.csv', index = False)