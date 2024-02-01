
import pandas as pd
import structman.lib.database.database as smdb


def score_consensus(df):
    """
    Rank templates based on majority consensus
    i.e. which has the most positions whose structural annotations agree with majority
    """
    # get frequency, count, unique, and top annotation
    cls = df.groupby(['Protein-ID', 'Position'])['Residue RIN simple classification'].describe()
    cls['consensus'] = cls['freq'] / cls['count']
    cls['majority_threshold'] = 1 / cls['unique']

    # ignore positions where unanimous or no majority (implicitly ignores positions with less than 3 template options)
    cls = cls.loc[cls['consensus'] > cls['majority_threshold']]

    # for each template, sum and rank the number of positions that match the consensus annotation
    df = df.set_index(['Protein-ID', 'Position']).join(cls['top'])
    df.dropna(subset=['top'], inplace=True)
    df['Score'] = df['top'] == df['Residue RIN simple classification']
    rank = df.groupby(['Protein-ID', 'PDB-ID', 'Chain'])['Score'].sum().to_frame()
    rank['Rank'] = rank.groupby(['Protein-ID']).transform('rank', ascending=False)
    rank.sort_values(by=['Protein-ID', 'Score'], ascending=False)

    return rank


def score_weighted_consensus(df):
    raise NotImplementedError('weighted consensus ranking method not implemented, please use unweighted consensus')


# called by structman_main
def suggest(session, config, full_annotations_file, outfile, ranking='consensus'):
    """
    Rank suggested templates based on ranking function
    """
    # map ranking method to corresponding function
    method = {
        'consensus': score_consensus,
        'weighted_consensus': score_weighted_consensus,
    }

    df = pd.read_csv(full_annotations_file, sep='\t', na_values=['-'])
    rank = method[ranking](df)
    rank = rank.reset_index()
    rank['Method'] = ranking

    # map names to db ids and insert ranked suggestions into database

    proteins = smdb.proteinsFromDb(session, config, with_alignments=True)
    protein_id_map = {name: proteins.get_protein_db_id(name) for name in proteins.get_protein_map()}
    structure_id_map = {name + chain: proteins.get_structure_db_id(name, chain) for name, chain in proteins.get_structures()}
    rank['Structure'] = (rank['PDB-ID'] + rank['Chain']).map(structure_id_map)  # structure db id
    rank['Protein'] = rank['Protein-ID'].map(protein_id_map)  # protein db id
    smdb.insert_suggestions(config, rank)  # store in database

    # save to file
    outcols = ['Protein-ID', 'PDB-ID', 'Chain', 'Method', 'Score', 'Rank']
    rank[outcols].sort_values(by=['Protein-ID', 'Method', 'Rank', 'PDB-ID']).to_csv(outfile, sep='\t', index=False)

    if config.verbosity >= 1:
        print(f"Templates ranked using `{ranking.replace('_', ' ')}`.\nSaved to {outfile}")

    return rank
