from setuptools import setup

setup(
  name = 'protclass',       
  packages = ['protclass'],   
  version = '1.0',      
  license='MIT',        
  description = 'Preprocessing and feature engineering for proteins and peptides',   
  author = 'Thomas Dorfer',                   
  author_email = 'thomas.a.dorfer@gmail.com',   
  url = 'https://github.com/thomasadorfer/ProtClass',   
  download_url = 'https://github.com/thomasadorfer/ProtClass/archive/v_01.tar.gz',  
  keywords = ['proteins', 'peptides', 'preprocessing', 'feature engineering', 'AA Index'], 
  install_requires=[            
          'validators',
          'beautifulsoup4',
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