<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/imgs/protlearn_logo.png" height="85" width="230">
</p>

<p align="center">
  A Python package for extracting protein sequence features
  <br>
  <a href="https://protlearn.readthedocs.io/en/latest/">Documentation</a>
  ·
  <a href="https://github.com/tadorfer/protlearn/issues/new?assignees=&labels=&template=feature_request.md&title=%5BNEW+FEATURE%5D">Request a feature</a>
  · 
  <a href="https://github.com/tadorfer/protlearn/issues/new?assignees=&labels=&template=bug_report.md&title=%5BBUG%5D">Report a bug</a>
  <br><br>
  <a href="https://travis-ci.org/tadorfer/protlearn"><img alt="Travis CI" src="https://img.shields.io/travis/tadorfer/protlearn"></a>
  <a href="https://codecov.io/gh/tadorfer/protlearn"><img alt="Codecov" src="https://codecov.io/gh/tadorfer/protlearn/branch/master/graph/badge.svg"></a>
  <a href="https://protlearn.readthedocs.io/en/latest/?badge=latest"><img alt="Docs" src="https://readthedocs.org/projects/protlearn/badge/?version=latest"></a> 
  <a href="https://pypi.org/project/protlearn/"><img alt="PyPI" src="https://img.shields.io/pypi/v/protlearn"></a>
  <a href="https://anaconda.org/conda-forge/protlearn"><img alt="Conda version" src="https://img.shields.io/conda/vn/conda-forge/protlearn.svg"></a>
  <a href="https://img.shields.io/pypi/pyversions/protlearn"><img alt="Python versions" src="https://img.shields.io/pypi/pyversions/protlearn"></a>  
  <a href="https://lbesson.mit-license.org/"><img alt="License" src="https://img.shields.io/badge/License-MIT-blue.svg"></a>   
</p>
<hr><br>

*protlearn* is a Python package for the feature extraction of amino acid sequences.
It is comprised of three stages - preprocessing, feature computation, and 
subsequent dimensionality reduction. Currently, the package is being maintained 
for Python versions 3.6, 3.7, and 3.8. 

## Overview

<p align="center">
  <img src="/imgs/protlearn_summary.png" height="430" width="624">
</p>

For more information on how to use it, please refer to the [documentation](https://protlearn.readthedocs.io/en/latest/).

## Installation

### Dependencies

- NumPy 
- Pandas 
- scikit-learn
- xgboost
- mlxtend
- biopython

### User Installation

#### PyPI

You can install _protlearn_ with `pip`:

```
$ pip install protlearn
```

#### Conda

You can install _protlearn_ with `conda`:

```
$ conda install -c conda-forge protlearn
```

## Authors

This package is maintained by [Thomas Dorfer](https://github.com/tadorfer).

## License

This package is licensed under the [MIT License](https://github.com/tadorfer/ProtLearn/blob/master/LICENSE).