import warnings
import numpy as np
import pandas as pd

from sklearn.preprocessing import MinMaxScaler

def scale_series(series,feature_range=(0, 1)):
    name = series.name
    scaler = MinMaxScaler(feature_range=feature_range)
    series_2d = series.values.reshape(-1, 1)
    scaled_series_2d = scaler.fit_transform(series_2d)
    scaled_series = pd.Series(scaled_series_2d.flatten(), index=series.index)
    scaled_series.name =name
    return scaled_series
    
def sort_by_terms_order(categorical_series, order_list, numerical_series):
    df = pd.DataFrame({
        'cat': categorical_series,
        'num': numerical_series
    })
    df = df[df['cat'].isin(order_list)]
    cat_type = pd.CategoricalDtype(categories=order_list, ordered=True)
    df['cat'] = df['cat'].astype(cat_type)
    sorted_df = df.sort_values(by=['cat', 'num'])
    return sorted_df.index
    
def df_fisher_chi2(clusters = pd.Series, response=pd.Series, R=False, NR=True):
    import pandas as pd
    from scipy.stats import fisher_exact, chi2_contingency
    from statsmodels.stats.multitest import multipletests
    df=pd.crosstab(clusters,response)
    df.insert(0,'Fisher_pv', 1)
    df.insert(1,'Chi2_pv', 1)

    for  i in df.index:
        nr, r = df[R].loc[i], df[NR].loc[i]
        nrj = df[R].sum() - nr
        rj = df[NR].sum() - r
        oddsratio, pvalue = fisher_exact([[nr, r],[nrj, rj]])  
        if pvalue > 1:
            pvalue = 1
        df.at[i,'Fisher_pv'] = pvalue

    for  i in df.index:
        nr, r = df[R].loc[i], df[NR].loc[i]
        nrj = df[R].sum() - nr
        rj = df[NR].sum() - r
        chi, pvalue, dof, exp = chi2_contingency([[r, nr],[rj, nrj]])  
        if pvalue > 1:
            pvalue = 1
        df.at[i, 'Chi2_pv'] = pvalue
    _, df['Fisher_pv'],_,_ = multipletests(df['Fisher_pv'],method='fdr_bh')
    _, df['Chi2_pv'],_,_ = multipletests(df['Chi2_pv'],method='fdr_bh')
    return df
    
class GeneSet(object):
    def __init__(self, name, descr, genes):
        self.name = name
        self.descr = descr
        self.genes = set(genes)
        self.genes_ordered = list(genes)

    def __str__(self):
        return '{}\t{}\t{}'.format(self.name, self.descr, '\t'.join(self.genes))


def read_gene_sets(gmt_file):
    """
    Return dict {geneset_name : GeneSet object}

    :param gmt_file: str, path to .gmt file
    :return: dict
    """
    gene_sets = {}
    with open(gmt_file) as handle:
        for line in handle:
            items = line.strip().split('\t')
            name = items[0].strip()
            description = items[1].strip()
            genes = set([gene.strip() for gene in items[2:]])
            gene_sets[name] = GeneSet(name, description, genes)

    return gene_sets


def ssgsea_score(ranks, genes):
    common_genes = list(set(genes).intersection(set(ranks.index)))
    if not len(common_genes):
        return pd.Series([0] * len(ranks.columns), index=ranks.columns)
    sranks = ranks.loc[common_genes]
    return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (len(ranks.index) - len(common_genes) + 1) / 2


def ssgsea_formula(data, gene_sets, rank_method='max'):
    """
    Return DataFrame with ssgsea scores
    Only overlapping genes will be analyzed

    :param data: pd.DataFrame, DataFrame with samples in columns and variables in rows
    :param gene_sets: dict, keys - processes, values - bioreactor.gsea.GeneSet
    :param rank_method: str, 'min' or 'max'.
    :return: pd.DataFrame, ssgsea scores, index - genesets, columns - patients
    """

    ranks = data.T.rank(method=rank_method, na_option='bottom')

    return pd.DataFrame({gs_name: ssgsea_score(ranks, gene_sets[gs_name].genes)
                         for gs_name in list(gene_sets.keys())})


def median_scale(data, clip=None):
    c_data = (data - data.median()) / data.mad()
    if clip is not None:
        return c_data.clip(-clip, clip)
    return c_data


def read_dataset(file, sep='\t', header=0, index_col=0, comment=None):
    return pd.read_csv(file, sep=sep, header=header, index_col=index_col,
                       na_values=['Na', 'NA', 'NAN'], comment=comment)


def item_series(item, indexed=None):
    """
    Creates a series filled with item with indexes from indexed (if Series-like) or numerical indexes (size=indexed)
    :param item: value for filling
    :param indexed:
    :return:
    """
    if indexed is not None:
        if hasattr(indexed, 'index'):
            return pd.Series([item] * len(indexed), index=indexed.index)
        elif type(indexed) is int and indexed > 0:
            return pd.Series([item] * indexed, index=np.arange(indexed))
    return pd.Series()


def to_common_samples(df_list=()):
    """
    Accepts a list of dataframes. Returns all dataframes with only intersecting indexes
    :param df_list: list of pd.DataFrame
    :return: pd.DataFrame
    """
    cs = set(df_list[0].index)
    for i in range(1, len(df_list)):
        cs = cs.intersection(df_list[i].index)

    if len(cs) < 1:
        warnings.warn('No common samples!')
    return [df_list[i].loc[list(cs)] for i in range(len(df_list))]


def cut_clustermap_tree(g, n_clusters=2, by_cols=True, name='Clusters'):
    """
    Cut clustermap into desired number of clusters. See scipy.cluster.hierarchy.cut_tree documentation.
    :param g:
    :param n_clusters:
    :param by_cols:
    :param name:
    :return: pd.Series
    """
    from scipy.cluster.hierarchy import cut_tree
    if by_cols:
        link = g.dendrogram_col.linkage
        index = g.data.columns
    else:
        link = g.dendrogram_row.linkage
        index = g.data.index

    return pd.Series(cut_tree(link, n_clusters=n_clusters)[:, 0], index=index, name=name) + 1


def pivot_vectors(vec1, vec2, na_label_1=None, na_label_2=None):
    """
    Aggregates 2 vectors into a table with amount of pairs (vec1.x, vec2.y) in a cell
    Both series must have same index.
    Else different indexes values will be counted in a_label_1/na_label_2 columns if specified or ignored
    :param vec1: pd.Series
    :param vec2: pd.Series
    :param na_label_1: How to name NA column
    :param na_label_2: How to name NA row
    :return: pivot table
    """

    name1 = str(vec1.name)
    if vec1.name is None:
        name1 = 'V1'

    name2 = str(vec2.name)
    if vec2.name is None:
        name2 = 'V2'

    if name1 == name2:
        name1 += '_1'
        name2 += '_2'

    sub_df = pd.DataFrame({name1: vec1,
                           name2: vec2})
    # FillNAs
    fill_dict = {}
    if na_label_1 is not None:
        fill_dict[name1] = na_label_1
    if na_label_2 is not None:
        fill_dict[name2] = na_label_2
    sub_df.fillna(value=fill_dict, inplace=True)

    sub_df = sub_df.assign(N=item_series(1, sub_df))

    return pd.pivot_table(data=sub_df, columns=name1,
                          index=name2, values='N', aggfunc=sum).fillna(0).astype(int)

def iterative_pca_outliers(df, return_labels=True):
    from sklearn.decomposition import PCA
    from scipy.stats import median_abs_deviation as mad
    pca = PCA(n_components=2)
    continue_search = True
    labels = []
    data = df.values
    index = df.index
    while continue_search:
        transformed_data = pca.fit_transform(data)
        medians = np.median(transformed_data, axis=0)
        mads = mad(transformed_data, scale=1/1.4826, axis=0)
        a = np.abs((transformed_data - medians) / mads)
        outliers_idx = np.any(a > 6, axis=1)
        if np.any(outliers_idx):
            for idx in index[outliers_idx]:
                print(f'{idx} is detected as an outlier')
                labels.append(idx)
            data = data[~outliers_idx]
            index = index[~outliers_idx]
        else:
            continue_search = False
            print('There are no outliers')
    if return_labels:
        return labels
