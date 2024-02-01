import pandas as pd
import structman.lib.database as smdb


def insert_suggestions(config, rank_df):
    insert_df = rank_df.copy()
    cols = ['Protein', 'Structure', 'Method', 'Score', 'Rank']
    check_dups = ['Protein', 'Structure', 'Method']  # combination of columns that should be unique

    # remove duplicates before inserting into db
    # first select all overlapping protein/structure ids already in database
    existing = smdb.select(config, cols, 'Suggestion', in_rows={col: insert_df[col].values for col in check_dups})
    existing = pd.DataFrame(existing, columns=cols)
    existing['dup'] = True  # add duplicate identifier column
    # merge to add duplicate identifier to new suggestions
    insert_df = insert_df.merge(existing, how='left', on=check_dups, suffixes=['', 'dup']).fillna({'dup': False})
    # drop duplicates before inserting
    insert_df = insert_df[~insert_df['dup']]

    smdb.insert('Suggestion', cols, insert_df[cols].values, config)  # insert non-duplicates into database
