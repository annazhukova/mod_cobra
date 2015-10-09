__author__ = 'anna'

from setuptools import setup, find_packages
from sys import version

if version < '2.2.3':
    from distutils.dist import DistributionMetadata

    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None


setup(name='mod_cobra',
      description='Combining several constraint-based analysis methods to analyse metabolic models.',
      long_description=open('README.md').read(),
      author='Anna Zhukova',
      author_email='zhutchok@gmail.com',
      version='0.1',
      packages=find_packages(),
      package_data={'mod_cobra.html': ['templates/*.html']},
      include_package_data=True,
      platform=['MacOS', 'Linux', 'Windows'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      download_url='https://github.com/annazhukova/mod_cobra',
      install_requires=['python-libsbml-experimental', 'numpy', 'cobra', 'mod_sbml', 'scipy', 'jinja2', 'igraph',
                        'louvain', 'sbml_vis']
      )

