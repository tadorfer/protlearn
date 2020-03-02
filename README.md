# Project description

[![PyPI](https://img.shields.io/pypi/v/ProtClass)](https://img.shields.io/pypi/v/ProtClass)
[![PyPI version fury.io](https://badge.fury.io/py/ansicolortags.svg)](https://pypi.org/project/protclass/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

## protclass

protclass is a Python module for preprocessing amino acid sequences (i.e. 
proteins and peptides) and subsequent feature engineering, both of which are
crucial steps to take prior to classification or regression problems.

The module is distributed under the MIT license and is being maintained by
Thomas Dorfer.

## Installation

### Dependencies

- Python 
- NumPy 
- Pandas 
- scikil-learn

### User Installation

```
$ pip install protclass
```

## Documentation

Currently, protclass is comprised of two preprocessing functions and four 
feature engineering functions.

### Preprocessing

- <i>txt_to_df</i>
- <i>integer_encoding</i>

#### txt_to_df

When working with protein or peptide sequences, the data almost always comes in
form of raw .txt files containing these sequences (typically one per row). 

This function converts these sequences into a Python-friendly Pandas DataFrame
with one column containing the sequences ['Sequence'] and a second column 
containing their corresponding class ['Label'] (this could correspond to the
protein family the sequence belongs, or whether a peptide is immunogenic or
not, etc.).

```
data_df = txt_to_df(raw_txt_file, integer_label)
```

#### integer_encoding

Machine learning algorithms can only handle numerical inputs. Therefore, the 
amino acid sequences need to be converted into numerical information, which is
achieved in the form of integers. 

This function converts amino acids of proteins or peptides into corresponding
integer values between 1-20. Zero, in this case, is reserved for padding these
sequences at the end to make them conform to a universal length (i.e. the 
length of the longest sequence in the dataset).

```
enc = integer_encode(data_df, padding=False)
```

### Feature engineering

- <i>length</i>
- <i>composition</i>
- <i>aaindex1</i>
- <i>aaindex3</i>

#### length

This function returns an n-dimensional array containing the lengths of all
sequences. The method can also be set to 'ohe', short for one-hot-encoding,
which leads to the generation of an array with n rows and the number of columns
corresponding to the longest sequence.

```
lengths = length(data_df, method='int')
```

#### composition

This function returns an array of shape (n, 20) containing the absolute or
relative frequencies of each amino acid that the sequence is comprised of.

```
comp = composition(data_df, method='relative', round_fraction=3)
```

#### aaindex1

This function computes the physicochemical properties of each amino acid 
comprising the sequence and returns the mean for each amino acid index.

```
aand1 = aaindex1(data_df, standardize='none')
```

#### aaindex3

This function computes the pairwise contact potentials between all amino acids
of a sequence and returns the mean value of each index.

``` 
aaind3 = aaindex3(data_df, standardize='none')
```