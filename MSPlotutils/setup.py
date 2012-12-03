from setuptools import setup, find_packages
setup(
    name='MSPlot',
    version="0.1.3.5j",
    description="Utilities for plotting 3D mass spec data, XIC etc from mzML data files",
    packages=['MSPlot'],
    install_requires=['matplotlib', 'pyteomics', 'pymzml'],
    author = "David Martin",
    author_email = "d.m.a.martin@dundee.ac.uk",
    )

