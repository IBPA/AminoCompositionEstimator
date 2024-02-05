# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='ace',
    version='1.0.0',
    description='A fast quantitative approach to estimate amino acid'
                'composition',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/IBPA/ACE',  # GitHub link.
    author='Fangzhou Li',
    author_email='fzli@ucdavis.edu',
    keywords=[
        'amino acid composition',
        'proteomics data',
        'quantitative analysis'
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'amino_acid_composition_estimate=ace.main:run']
    },
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.19.3',
        'pandas>=1.1.5',
        'notebook>=6.1.5'
    ]
)
