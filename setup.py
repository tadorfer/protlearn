from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'protlearn',       
  packages = find_packages(exclude=["tests.*", "tests"]), 
  package_data={'protlearn': ['features/data/*.csv']},  
  version = '0.0.2',      
  license='MIT',        
  description = 'A Python package for extracting protein sequence features', 
  long_description_content_type='text/markdown',
  long_description=long_description,  
  author = 'Thomas Dorfer',                   
  author_email = 'thomas.a.dorfer@gmail.com',   
  url = 'https://github.com/tadorfer/protlearn',   
  download_url = 'https://github.com/tadorfer/protlearn/archive/v0.0.2.tar.gz',  
  keywords = ['amino acids', 'proteins', 'peptides', 'preprocessing', 'feature engineering', 'dimensionality reduction', 'machine learning'], 
  setup_requires = ['wheel'],
  install_requires=[            
          'numpy',
          'pandas',
          'scikit-learn',
          'xgboost',
          'mlxtend',
          'biopython'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      
    'Intended Audience :: Science/Research',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ]
)
