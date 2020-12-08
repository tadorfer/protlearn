.. _dimensionality_reduction:
.. |br| raw:: html

   <br />

Dimensionality Reduction & Feature Selection
============================================

univariate_filter
-----------------

.. code-block:: text

    protlearn.dimreduction.univariate_filter(X, y, *, method='f_test', top=10)

Univariate feature selection.

This function returns the features selected by univariate filtering after 
examining each feature individually and determining the strength of its 
relationship with the response variable. Here, three statistical tests can 
be chosen: f-test, chi-squared, and mutual information.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre)
    Feature matrix.

y: ndarray of shape (n_samples,)
    Response variables.

method: string, default='f_test'
    'f_test' : ANOVA f-scores |br|
    'chi2' : Chi-squared statistics |br|
    'mutual_info' : Mutual information

top: int, default=10
    Number of top features to select.
    
Returns
#######

arr: ndarray of shape (n_samples, top)
    Array containing the top features.

Examples
########

.. code-block:: python

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
    
correlation
-----------

.. code-block:: text

    protlearn.dimreduction.correlation(X, thres=.9)

Pearson correlation.

This function returns the features whose Pearson correlation with one 
another is below a specified threshold, thus circumventing the problem of 
multicollinearity.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix.

thres: float, default=.9
    Features whose correlation coefficient is higher than this threshold 
    value will be removed.

Returns
#######

arr:  ndarray of shape (n_samples, n_features_post)
    Array containing features that correlate below the threshold with one 
    another.

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import correlation
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = correlation(features, thres=.99)
    >>> reduced.shape
    (3, 12)

lasso
-----

.. code-block:: text

    protlearn.dimreduction.lasso(X, y, C=1.0)

Lasso (L1) regularization.

Linear Model trained with L1 prior as regularizer. 

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre)
    Feature matrix.

y: ndarray of shape (n_samples,)
    Response variables.

C: float, default=1.0
    Inverse of regularization strength.
    
Returns
#######

arr : ndarray of shape (n_samples, n_features_post)
    Array containing lasso-reduced features.

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import lasso
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = lasso(features, labels)
    >>> reduced.shape
    (3, 2)

tree_importance
---------------

.. code-block:: text

    protlearn.dimreduction.tree_importance(X, y, *, clf=None, method='random_forest', top=None, n_estimators=100, max_depth=None, importance_type='gain')

Decision tree feature importance.

This function returns the features that were selected as important by 
decision tree algorithms such as Random Forest and XGBoost. 

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix.

y: ndarray of shape (n_samples,)
    Response variables.

clf: object or None, default=None
    Customized classifier.

method: string, default='random_forest'
    'random_forest' : Random Forest Classifier
    'xgboost' : XGBoost Classifier

top: int or None, default=None
    Number of top features to select.

n_iterations: int, default=3
    Number of iterations.
    
n_estimators: int or None, default=2
    Number of trees in the forest.
    
max_depth: int or None, default=None
    Maximum depth of the tree.
    
importance_type: string, default='gain'
    For XGBoost only: |br|
    'gain' : average gain of splits which use the feature |br|
    'weight' : number of times the a feature appears in the tree |br|
    'cover' : average coverage of splits which use the feature |br|
    'total_gain' : Total gain |br|
    'total_cover' : Total cover

Returns
#######

arr:  ndarray of shape (n_samples, top)
    Array containing the top features based on tree-importance.

indices:  ndarray
    Indices indicating the position of the selected feature in the input vector.


Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import tree_importance
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced, indices = tree_importance(features, labels, top=10)
    >>> reduced.shape
    (3, 10)
    >>> indices
    array([249, 514, 4, 155, 182,  82, 214, 405, 140, 364])

sequential
----------

.. code-block:: text

    protlearn.dimreduction.sequential(X, y, *, estimator, direction='forward', n_features=10, cv=0)

Sequential feature selection.

Sequential feature selection algorithms are a family of greedy search 
algorithms that are used to reduce an initial d-dimensional feature space 
to a k-dimensional feature subspace where k < d. These algorithms remove or 
add one feature at a time based on the classifier performance until a 
feature subset of the desired size k is reached.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix.

y: labels, ndarray of shape (n_samples,)
    Response variables.

estimator: object
    Classifier - must include \coef_ or \feature_importances_ attribute.
    
direction: string, default='forward'
    Direction of sequential model, can be 'forward' or 'backward'.

n_features: int, default=None
    Number of features to select.
    
cv: int, default=0
    Number of cross-validation steps.

Returns
#######

arr:  ndarray of shape (n_samples, n_features)
    Array containing features selected by the sequential models.

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import sequential
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> rf = RandomForestClassifier()
    >>> reduced = sequential(features, labels, rf, n_features=10)
    >>> reduced.shape
    (3, 10)

rfe
---

.. code-block:: text

    protlearn.dimreduction.rfe(X, y, *, estimator, n_features=None, step=1)

Recursive feature elimination.

This function selects features by recursively considering smaller and 
smaller feature subsets. First, the estimator is trained on the initial 
feature matrix and the importance of each feature is obtained through a 
\coef_ or a \feature_importances_ attribute. Subsequently, the least 
important features are pruned from the current feature subset. This is 
repeated recursively on the pruned subset until the desired number of 
features is eventually reached.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix.

y: labels, ndarray of shape (n_samples,)
    Response variables.

estimator: object
    Classifier - must include \coef_ or \feature_importances_ attribute.
    
n_features: int or None, default=None
    Number of features to select. If ``None``, half of the features are selected.
    
step: int, default=1
    Number of features to remove at each iteration.

Returns
#######

arr:  ndarray of shape (n_samples, n_features)
    Array containing the RFE-selected features.

ranking: ndarray of shape (n_features_pre,)
    Ranking of the features (with 1 being the best).

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import rfe
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> rf = RandomForestClassifier()
    >>> reduced, _ = rfe(features, labels, rf, n_features=10, step=5)
    >>> reduced.shape
    (3, 10)

pca
---

.. code-block:: text

    protlearn.dimreduction.pca(X, *, thres=.9, whiten=False)

Principal component analysis.

PCA is defined as an orthogonal linear transformation that transforms the 
data to a new coordinate system such that the greatest variance by some 
scalar projection of the data comes to lie on the first coordinate (called 
the first principal component), the second greatest variance on the second 
coordinate, and so on.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre)
    Feature matrix. 

thres: float, default=.9
    Specify the desired explained variance.

Returns
#######

arr:  ndarray of shape (n_samples, n_features_post)
    Array containing the PCA components comprising the specified variance.

Notes
#####

For the output to be meaningful, the number of samples should be larger than
the number of features.

Examples
########

.. code-block:: python

    >>> from protlearn.dimreduction import pca
    >>> features.shape #from a larger dataset (not shown here)
    (1000, 575)
    >>> reduced = pca(features, thres=.9)
    (1000, 32)

lda 
---

.. code-block:: text

    protlearn.dimreduction.lda(X, y, *, solver='svd', shrinkage=None, n_components=None)

Linear discriminant analysis.

This function reduces the dimensionality of the input by projecting it to 
the most discriminative directions.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix. 

y: ndarray of shape (n_samples,)
    Response variables.

solver: string, default='svd'
    'svd' : Singular value decomposition |br|
    'lsqr' : Least squares solution |br|
    'eigen' : Eigenvalue decomposition
    
shrinkage: string, float, or None, default=None
    Shrinkage parameter. |br|
    None : no shrinkage |br|
    'auto' : automatic shrinkage using the Ledoit-Wolf lemma |br|
    float between 0 and 1: fixed shrinkage parameter
    
n_components: int or None, default=None
    Number of components for dimensionality reduction. This parameter 
    cannot be larger than min(n_features, n_classes - 1).

Returns
#######

arr:  ndarray of shape (n_samples, n_features_post)
    Array containing the LDA-transformed features.

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import lda
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> labels = [1., 0., 0.]
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = lda(features, labels, n_components=1)
    >>> reduced.shape
    (3, 1)

tsne 
----

.. code-block:: text

    protlearn.dimreduction.tsne(X, *, n_components=2, perplexity=30, prior_pca=True, pca_components=50)

t-distributed stochastic neighbor embedding.

t-SNE converts similarities between data points to joint probabilities and 
tries to minimize the Kullback-Leibler divergence between the joint 
probabilities of the low-dimensional embedding and the high-dimensional data.

Parameters
##########

X: ndarray of shape (n_samples, n_features_pre) 
    Feature matrix.
    
n_components: int or None, default=2
    Dimension of embedded space.
    
perplexity: int, default=30
    Related to the number of nearest neighbors that is used in other 
    manifold learning algorithms. Should be between 5 and 50. Larger 
    datasets require larger perplexity.
    
prior_pca: bool, default=True
    It is recommended to reduce dimensionality before running t-SNE to 
    decrease computation time and noise.
    
pca_components: int, default=50
    Dimension of PCA-preprocessed data that will serve as input to t-SNE.

Returns
#######

arr:  ndarray of shape (n_samples, n_components)
    Array containing the t-SNE-transformed features.

Examples
########

.. code-block:: python

    >>> import numpy as np
    >>> from protlearn.features import aac, aaindex1, ngram
    >>> from protlearn.dimreduction import tsne
    >>> seqs = ['ARKLY', 'EERKPGL', 'PGPGEERNLY']
    >>> comp, _ = aac(seqs)
    >>> aaind, _ = aaindex1(seqs)
    >>> ng, _ = ngram(seqs)
    >>> features = np.concatenate([comp, aaind, ng], axis=1)
    >>> features.shape
    (3, 575)
    >>> reduced = tsne(features, pca_components=3)
    >>> reduced.shape
    (3, 2)
