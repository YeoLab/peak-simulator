from distutils.core import Extension
from setuptools import setup, find_packages

peaks = Extension("src/peaks", sources = ['src/peaksmodule.cc'])
setup(
    name = "FindPeaks",
    version = "0.1",
    packages = find_packages(),
    ext_modules = [peaks],
    #Project Dependences go here

    #data and scripts go here

    #metadata for upload to PyPI
    author = "Michael Lovci and Gabriel Pratt",
    author_email = "lovci@ucsd.edu",
    description = "A set of scripts for calling peaks on CLIP-seq data",
    license = "TBD",
    keywords = "CLIP-seq, peaks, bioinformatics",
    url = "http://yeolab.ucsd.edu/",
    
    #Other stuff I feel like including here
    include_package_data = True,
    #zip_safe = True #True I think
)
