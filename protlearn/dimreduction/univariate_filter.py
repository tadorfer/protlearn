# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif, chi2, mutual_info_classif
from sklearn.preprocessing import MinMaxScaler

def univariate_filter(X, y, *, method='f_test', top=10):
    """Univariate feature selection.

    This function returns the features selected by univariate filtering after 
    examining each feature individually and determining the strength of its 
    relationship with the response variable. Here, three statistical tests can 
    be chosen: f-test, chi-squared, and mutual information.
    
    Parameters
    ----------
    
    X : ndarray of shape (n_samples, n_features_pre)
        Feature matrix.
    
    y : ndarray of shape (n_samples,)
        Response variables.
    
    method : string, default='f_test'
        'f_test' : ANOVA f-scores
        'chi2' : Chi-squared statistics
        'mutual_info' : Mutual information
    
    top : int, default=10
        Number of top features to select.
        
    Returns
    -------
    
    arr : ndarray of shape (n_samples, top)
        Array containing the top features.
    
    Examples
    --------

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import univariate_filter
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = univariate_filter(features, labels, method='f_test', top=10)
    >>> reduced.shape
    (3, 10)
    
    """  
    
    if method == 'f_test':
        arr = SelectKBest(f_classif, k=top).fit_transform(X, y)
    elif method == 'chi2':
        # only non-negative features
        try:
            arr = SelectKBest(chi2, k=top).fit_transform(X, y)
        except:
            X = MinMaxScaler().fit_transform(X)
            arr = SelectKBest(chi2, k=top).fit_transform(X, y)
    elif method == 'mutual_info':
        arr = SelectKBest(mutual_info_classif, k=top).fit_transform(X, y)
    
    return arr