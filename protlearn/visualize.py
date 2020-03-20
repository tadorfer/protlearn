# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import pandas as pd
import seaborn as sns
import matplotlib as mpl
from collections import Counter
from feature_engineering import ngram_composition


# set sns style settings
sns.set(font_scale=1.3)
sns.set_style('white', {'axes.spines.top': False, 'axes.spines.right': False})

def viz_length(X, method='absolute', sort=True, get_data=False, plot=True):
    """Bar plot of the length of a sequence or set of sequences. 

    While it works on single sequences, this function is really only useful
    when investigating a set of sequences to show how many sequences are of
    which length.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    method : string, default='absolute'

        'absolute' : visualize absolute frequency of sequence length
        'relative' : visualize relative frequency of sequence length

    sort : bool, default=True
        Sort by descending lengths.

    get_data : bool, default=True
        Return the corresponding dataframe.

    plot : bool, default=True
        Show bar plot.

    Returns
    -------

    plot: seaborn barplot showing the frequency of sequence lengths

    df_lengths : Pandas DataFrame of shape (1, n_unique_lengths) if 
                 get_data = True.
        The column names are the lengths of the sequences. 

    Notes
    -----
    For comparing the sequence lengths of various datasets, it is recommended to
    pass method='relative', which returns the frequency of each length divided 
    by the total number of sequences. 

    """

    len_all = [len(seq) for seq in X['Sequence']]
    len_unique = list(set(len_all))
    len_dict = {str(len_unique[i]): len_all.count(x) for i, x in\
                enumerate(len_unique)}
    lengths = list(len_dict.values())
    
    columns = [i for i in list(len_dict.keys())]
    df = pd.DataFrame(data=[lengths], columns=columns, index=['Frequency'])
    df = df.reindex(sorted(df.columns), axis=1)
    
    if method == 'relative':
        df.iloc[0,:] /= (df.sum(axis=1)/100)[0]
    
    if sort == True:
        df = df.sort_values('Frequency', axis=1, ascending=False).\
                            reset_index(drop=True)
        df.index = ['Frequency']
    
    if plot == True:
        # plotting
        if method == 'relative':
            ylabel = 'Relative Frequency [%]'
        elif method == 'absolute':
            ylabel = 'Absolute Frequency'

        ax = sns.barplot(x=df.columns, y=df.iloc[0], 
                         palette=mpl.cm.ScalarMappable(cmap='coolwarm').\
                         to_rgba(df.iloc[0]), order=df.columns)
        ax.set_xlabel('Length [amino acids]', fontsize=15, labelpad=10)
        ax.set_ylabel(ylabel, fontsize=15, labelpad=10)
    
    if get_data == True:
        return df


def viz_composition(X, method='absolute', sort=True, get_data=False, plot=True):
    """Bar plot of the amino acid composition of a sequence or set of sequences.

    This function can be applied to single sequences or a dataset comprised of
    multiple sequences, in which case it will show the amino acid composition 
    of the entire dataset as a whole.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    method : string, default='absolute'

        'absolute' : visualize absolute frequency of amino acid composition
        'relative' : visualize relative frequency of amino acid composition

    sort : bool, default=True
        Sort by descending frequency.

    get_data : bool, default=True
        Return the corresponding dataframe.

    plot : bool, default=True
        Show bar plot.

    Returns
    -------

    plot: seaborn barplot showing amino acid composition

    df_comp : Pandas DataFrame of shape (20, 2) if get_data = True
        The two columns contain the amino acids and their frequency. 

    Notes
    -----
    For comparing the composition of various datasets, it is recommended to pass
    method='relative', which returns the frequency of each amino acid divided by
    the total length of all sequences. 

    """

    all_aa = [aa for seq in X['Sequence'] for aa in seq]
    aa_dict= Counter(all_aa)
    df = pd.DataFrame({'Amino Acid': list(aa_dict.keys()),\
                       'Frequency': list(aa_dict.values())})
    df = df.sort_values(by='Amino Acid').reset_index(drop=True)

    if sort == True:
        df = df.sort_values(by='Frequency',\
                            ascending=False).reset_index(drop=True)
    
    if method=='relative':
        df['Frequency'] /= sum(df['Frequency'])/100
        
    if plot==True:
        # plotting
        ax = sns.barplot(x='Amino Acid', 
                         y='Frequency', 
                         data=df, 
                         palette=mpl.cm.ScalarMappable(cmap='coolwarm').\
                         to_rgba(df['Frequency']),
                         order=df['Amino Acid'])

        ax.set_xlabel('Amino Acid',
                      fontsize=15, 
                      labelpad = 10)
        
        if method=='relative':
            ylabel = 'Relative Frequency [%]'
        elif method=='absolute':
            ylabel = 'Absolute Frequency'
            
        ax.set_ylabel(ylabel, fontsize=15, labelpad = 10)
        
    if get_data == True:
        return df


def viz_ngram(X, ngram=2, top=100, method='absolute', sort=True, get_data=False, 
              plot=True, xtick_rotation=0):
    """Bar plot of ngram composition.

    This function can be applied to single sequences or a dataset comprised of
    multiple sequences. In the latter case, it will show the ngram composition
    of all sequences in the dataset as a whole.

    Parameters
    ----------

    X : Pandas DataFrame 
        The column containing protein or peptide sequences must be labeled
        'Sequence'.

    ngram : int, default=2
        Integer denoting the desired n-gram composition.
        
        2 : dipeptide composition
        3 : tripepitde composition
        4 : quadpeptide composition
        
    top : int, default=100
        Display top X percent.

    method : string, default='absolute'

        'absolute' : visualize absolute ngram composition
        'relative' : visualize relative ngram composition

    sort : bool, default=True
        Sort by descending frequency.

    get_data : bool, default=True
        Return the corresponding dataframe.

    plot : bool, default=True
        Show bar plot.

    xtick_rotation : int, default=0
        If the plot is comprised of many ngram combination, the xticklabels can
        be rotated to not overlap with their neighbors.

    Returns
    -------

    plot: seaborn barplot showing ngram composition

    df_ngram : Pandas DataFrame of shape (1, unique_ngram_compositions) if
               get_data = True
        The column names are the ngram combinations. 

    Notes
    -----
    For comparing the ngram composition of various datasets, it is recommended 
    to pass method='relative', which returns the frequency of each ngram
    combination divided by the total ngram combinations. 

    """
    # compute ngram composition and get sum
    ng = ngram_composition(X, ngram)
    ng_sum = ng.sum()
    
    # build dataframe
    df = pd.DataFrame(data=[ng_sum], columns=ng_sum.index, index=['Frequency'])
    
    if method == 'relative':
        df.iloc[0,:] /= (df.sum(axis=1)/100)[0]
        
    if sort == True:
        df = df.sort_values('Frequency', axis=1, ascending=False).\
                            reset_index(drop=True)
        df.index = ['Frequency']

    if plot == True:
        # plotting
        if method == 'relative':
            ylabel = 'Relative Frequency [%]'
        elif method == 'absolute':
            ylabel = 'Absolute Frequency'

        if ngram == 2:
            xlabel = 'Dipeptide Composition'
        elif ngram == 3:
            xlabel = 'Tripeptide Composition'
        elif ngram ==4:
            xlabel = 'Quadpeptide Composition'
            
        first = round(len(df.columns) * (top/100))
        
        ax = sns.barplot(x=df.columns[:first],
                         y=df.iloc[0, :first],
                         palette=mpl.cm.ScalarMappable(cmap='coolwarm').\
                         to_rgba(df.iloc[0]), order=df.columns[:first])
        ax.set_xlabel(xlabel, fontsize=15, labelpad = 10)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=xtick_rotation)
        ax.set_ylabel(ylabel, fontsize=15, labelpad = 10)
    
    if get_data == True:
        return df