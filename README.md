[![Travis (.org)](https://img.shields.io/travis/tadorfer/ProtLearn)](https://travis-ci.org/tadorfer/ProtLearn)
[![PyPI](https://img.shields.io/pypi/v/ProtLearn)](https://pypi.org/project/protlearn/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/ProtLearn)](https://img.shields.io/pypi/pyversions/ProtLearn)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/protlearn)](https://pypistats.org/packages/protlearn)

# protlearn

protlearn is a Python package for preprocessing amino acid sequences (i.e. 
proteins and peptides) and subsequent feature engineering. Its functions are
particularly suited for the preparation of classification and regression tasks.

## Installation

### Dependencies

- NumPy 
- Pandas 
- scikit-learn
- seaborn
- matplotlib

### User Installation

```
$ pip install protlearn
```

## Documentation

* [Preprocessing](#preprocessing)
    - [txt_to_df](#txt_to_df)
    - [integer_encode](#integer_encode)
* [Feature Engineering](#feature-engineering)
    - [length](#length)
    - [composition](#composition)
    - [aaindex1](#aaindex1)
    - [aaindex2](#aaindex2)
    - [aaindex3](#aaindex3)
    - [ngram_composition](#ngram_composition)
* [Visualization](#visualization)  
    - [viz_length](#viz_length)
    - [viz_composition](#viz_composition)
    - [viz_ngram](viz_ngram)


### Preprocessing

#### `txt_to_df`

This function converts sequences from a raw `.txt file` into a Python-friendly 
Pandas DataFrame with the column name 'Sequence'. For classification tasks, it
provides the option of including an additional class column with the name 
'Label' by passing an integer value (denoting the class) to the function 
argument `label`. 

<br>

<b>Example:</b>

```python
from protlearn import txt_to_df

df = txt_to_df(test_seq.txt)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/text_to_df.png" height="260" width="460">
</p>

For more information --> `help(txt_to_df`)

<br>

#### `integer_encode`

This function converts amino acid sequences into corresponding integer values 
between 1-20. Zero, in this case, is reserved for the optional `padding` of 
these sequences at the end to make them conform to a universal length (i.e. the 
length of the longest sequence in the dataset).

The following amino acid order is used for conversion (1-20):

A C D E F G H I K L M N P Q R S T V W Y

<br>

<b>Example:</b>

```python
from protlearn import txt_to_df, integer_encode

df = txt_to_df(test_seq.txt)
enc = integer_encode(df, padding=True)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/integer_encode.png" height="220" width="560">
</p>

If `padding=True`, a numpy array of shape (n_samples, longest_sequence) will be
returned. Otherwise, a numpy array of shape (n_samples, ) containing each 
integer-encoded sequence as separate numpy arrays will be returned.

For more information --> `help(integer_encode)`

<br>

### Feature engineering

#### `length`

This function returns an n-dimensional array containing the length of all
sequences. The `method` can also be set to 'ohe', short for one-hot-encoding,
which leads to the generation of an array with n rows and the number of columns
corresponding to the number of unique sequence lengths.

<b>Example:</b>

```python
from protlearn import txt_to_df, length

df = txt_to_df(test_seq.txt)
lengths = length(df)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/length.png" height="300" width="460">
</p>

This illustration shows that, if `method='ohe'`, the columns correspond to the 
unique lengths of the sequences (in order). In this case, there is no sequence 
with length 8, so the columns correspond to sequence lengths 6, 7, and 9. 

For more information --> `help(length)`

<br>

#### `composition`

This function returns an array of shape (n_samples, n_unique_amino_acids) 
containing the absolute frequencies of each amino acid that the sequence is 
comprised of. If `method='relative'`, the absolute count of each amino acid is
divided by the sequence length and returned as a fraction, whose number of 
decimals can be chosen with the argument `round_fraction`. 

<b>Example:</b>

```python
from protlearn import txt_to_df, composition

df = txt_to_df(test_seq.txt)
comp = composition(df, method='absolute')
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/composition.png" height="250" width="590">
</p>

This illustration shows the absolute frequency of amino acids of each input
sequence. If a particular amino acid is not present in all of the input 
sequences, its column will not be returned to avoid all-zero columns. Therefore,
the number of columns of the returned dataframe is not always 20, but can vary.

For more information --> `help(composition)`

<br>

#### `aaindex1`

This function computes the physicochemical properties of each amino acid in the
sequence and returns the mean of each index per sequence. Currently, ver.9.2 
(release Feb, 2017) contains 566 indices. However, due to 13 of these indices
containing NaNs, the returned dataframe will have a column size of 553-566, 
depending on the number of amino acids in the sequence whose index is NaN, in
which case this index column will be removed. More information on the AAindex1 
can be found on [GenomeNet Database Resources](https://www.genome.jp/aaindex/).

<b>Example:</b>

```python
from protlearn import txt_to_df, aaindex1

df = txt_to_df(test_seq.txt)
aand1 = aaindex1(df, standardize='none')
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/aaindex1.png" height="250" width="760">
</p>

This illustration shows how AAIndex1 is computed, using 'ARN' as a sample 
sequence. It is highly recommended to `standardize` the data if the resulting
dataframe is intended to serve as input to a classifier/regressor.

For more information --> `help(aaindex1)`

<br>

#### `aaindex2`

This function computes the substitution matrices of each amino acid pair of a 
sequence and returns the mean of all substitution scores per sequence. 
Currently, ver.9.2 (release Feb, 2017) contains 94 indices. However, due to two
of these indices containing NaNs, the column size of the returned dataframe can
vary between 92-94, depending on whether any of the sequences contains an amino
acid pair whose index is NaN, in which case the entire column will be removed.
Furthermore, some of the square and rectangular matrices originally include gaps 
and various metrics for cysteines (disulfide-bonded and free), which are removed
for all computations (only the standard 20 amino acids are used here). More 
information on the AAindex2 can be found on [GenomeNet Database Resources](https://www.genome.jp/aaindex/).

<b>Example:</b>

```python
from protlearn import txt_to_df, aaindex2

df = txt_to_df(test_seq.txt)
aaind2 = aaindex2(df, standardize='none')
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/aaindex2.png" height="420" width="760">
</p>

This illustration shows how AAIndex2 is computed, using 'ARN' as a sample 
sequence. It also shows that this index contains both lower triangular as well 
as square matrices. It is highly recommended to `standardize` the data if the
resulting dataframe is intended to serve as input to a classifier/regressor.

For more information --> `help(aaindex2)`

<br>

#### `aaindex3`

This function computes the pairwise contact potentials between all amino acid
pairs in a sequence and returns the mean of all contact potentials per sequence.
Currently, ver.9.2 (release Feb, 2017) contains 47 indices. However, due to four
of these indices containing NaNs, the column size of the returned dataframe can
vary between 43-47, depending on whether any of the sequences contains an amino
acid pair whose index is NaN, in which case the entire column will be removed.
More information on the AAindex3 can be found on [GenomeNet Database Resources](https://www.genome.jp/aaindex/).

<b>Example:</b>

```python
from protlearn import txt_to_df, aaindex3

df = txt_to_df(test_seq.txt)
aaind3 = aaindex3(df, standardize='none')
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/aaindex3.png" height="420" width="760">
</p>

This illustration shows how AAIndex3 is computed, using 'ARN' as a sample 
sequence. As with AAIndex2, is also contains both lower triangular as well as
square matrices. It is highly recommended to `standardize` the data if the
resulting dataframe is intended to serve as input to a classifier/regressor.

For more information --> `help(aaindex3)`

<br>

#### `ngram_composition`

This function computes the di-, tri-, or quadpeptide composition of any given
amino acid sequence. 

<b>Example:</b>

```python
from protlearn import txt_to_df, ngram_composition

df = txt_to_df(test_seq.txt)
ngram = ngram_composition(df)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/ngram.png" height="350" width="450">
</p>

This illustration shows the di-, tri-, and quadpeptide composition of the 
sequence 'ALLLFTY'. As can be observed, not all possible combinations (20^ngram)
are returned, but only those whose values are non-zero.

For more information --> `help(ngram_composition)`

<br>

### Visualization

The following bar plots use a `coolwarm` color palette, meaning that there is a
gradient from higher bars (red) to lower bars (blue).

#### `viz_length`

This function returns a bar plot displaying the frequency of sequence lengths
in descending order.

<b>Example:</b>

```python
from protlearn import txt_to_df, viz_length

df = txt_to_df(test_seq.txt)
viz_length(df)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/viz_length.png" height="240" width="570">
</p>

This illustration shows the frequency of all lengths of sequences found in the 
`test_seq.txt` file. The data used for this plot can also be retrieved by 
passing the argument `get_data=True`. In this case, the function should be
called as follows: `lengths = viz_length(df, get_data=True)`.

For more information --> `help(viz_length)`

<br>

#### `viz_composition`

This function returns a bar plot displaying the frequency of amino acids
present in the dataset (in descending order).

<b>Example:</b>

```python
from protlearn import txt_to_df, viz_composition

df = txt_to_df(test_seq.txt)
viz_composition(df)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/viz_composition.png" height="250" width="600">
</p>

This illustration shows the frequency of all amino acids of sequences found in 
the `test_seq.txt` file. The data used for this plot can also be retrieved by 
passing the argument `get_data=True`. In this case, the function should be
called as follows: `comp = viz_composition(df, get_data=True)`.

For more information --> `help(viz_composition)`

<br>

#### `viz_ngram`

This function returns a bar plot displaying the frequency of the selected ngram
composition (di-, tri-, or quadpeptide) of the sequence (in descending order).

<b>Example:</b>

```python
from protlearn import txt_to_df, viz_ngram

df = txt_to_df(test_seq.txt)
viz_ngram(df, ngram=2, top=20)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/dems/viz_ngram.png" height="250" width="620">
</p>

This illustration shows the frequency of the top 20% (`top=20`) of dipeptide 
combinations (`ngram=2`) of sequences found in the `test_seq.txt` file. The 
data used for this plot can also be retrieved by passing the argument 
`get_data=True`. In this case, the function should be called as follows: 
`ng = viz_ngram(df, get_data=True)`.

For more information --> `help(viz_ngram)`

<br>

## Authors

This package is maintained by [Thomas Dorfer](https://github.com/tadorfer)

## License

This package is licensed under the [MIT License](https://github.com/tadorfer/ProtLearn/blob/master/LICENSE).