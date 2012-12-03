from setuptools import setup, find_packages
setup(
    name='Unimod',
    version="0.1.0",
    description="simple wrapper around the Unimod amino acid modificatins database",
    packages=find_packages(),
    install_requires=[''],
    package_data= {'':['unimod.xml'],},
    author = "David Martin",
    author_email = "d.m.a.martin@dundee.ac.uk",
    )

