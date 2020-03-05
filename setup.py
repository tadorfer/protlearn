from setuptools import setup
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'protclass',       
  packages = ['protclass'],   
  version = '1.4',      
  license='MIT',        
  description = 'Preprocessing and feature engineering for proteins and peptides prior to classification', 
  long_description_content_type='text/markdown',
  long_description=long_description,  
  author = 'Thomas Dorfer',                   
  author_email = 'thomas.a.dorfer@gmail.com',   
  url = 'https://github.com/tadorfer/ProtClass',   
  download_url = 'https://github.com/tadorfer/ProtClass/archive/v1.4.tar.gz',  
  keywords = ['proteins', 'peptides', 'preprocessing', 'feature engineering', 'AA Index'], 
  setup_requires = ['wheel'],
  include_package_data=True,
  install_requires=[            
          'numpy',
          'pandas',
          'cython',
          'scikit-learn'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Science/Research',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ],
)