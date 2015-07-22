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
      author_email='anna.zhukova@ibgc.cnrs.fr',
      version='0.1',
      packages=find_packages(),
      package_data={'mod_cobra.gibbs': ['data/*.csv']},
      include_package_data=True,
      license='LICENSE',
      platform=['MacOS', 'Linux', 'Windows'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: CeCILL',
          'Topic :: Systems Biology',
          'Topic :: Software Development',
      ],
      download_url='https://github.com/annazhukova/mod_cobra',
      install_requires=['openpyxl', 'python-libsbml-experimental', 'numpy', 'matplotlib', 'gmpy', 'cobra', 'networkx',
                        'mod_sbml', 'tulip', 'geojson', 'sympy', 'jinja2', 'libsbgnpy', 'tarjan']
)
