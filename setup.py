import os
from setuptools import setup, find_packages
import numpy
from pprocess import __version__

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='pprocess',
    version=__version__,
    description='A tool to automate resource descriptors of catalogues',
    url='https://github.com/epugaant/pprocess',
    author='E. Puga',
    author_email='elena@pugaantolin.eu',
    license='MIT',
    packages=find_packages(),
    long_description=read('README.md'),
    include_dirs=[numpy.get_include()],
    include_package_data = True,
    zip_safe=False,
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.x',
        ],
    keywords='Legacy Archive Astronomy Space Mission',
    install_requires=read('requirements.txt')
    )