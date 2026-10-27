"""
PyFLUPS - Python API for FLUPS

Setup script for installation
"""

from setuptools import setup, find_packages
import os

# Read README for long description
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

setup(
    name='pyflups',
    version='1.0.0',
    author='Thomas Gillis, Denis-Gabriel Caprace, and contributors',
    author_email='',
    description='Python API for FLUPS - Fourier-based Library of Unbounded Poisson Solvers',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/vortexlab-uclouvain/flups',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.18.0',
        'mpi4py>=3.0.0',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-mpi',
        ],
        'docs': [
            'sphinx>=4.0',
            'sphinx-rtd-theme',
        ],
    },
    project_urls={
        'Bug Reports': 'https://github.com/vortexlab-uclouvain/flups/issues',
        'Source': 'https://github.com/vortexlab-uclouvain/flups',
    },
)
