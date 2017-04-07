# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('yaps2/version.py') as f:
    exec(f.read())

setup(
    name='yaps2',
    version=__version__,
    description='Yet Another Pipeline Setup',
    long_description=readme,
    author='Indraniel Das',
    author_email='idas@wustl.edu',
    license=license,
    url='https://github.com/indraniel/yaps2',
    dependency_links=[
        'https://github.com/indraniel/COSMOS2/tarball/enable-lsf-rebase#egg=cosmos-wfm-2.0.10.indraniel4'
        'https://github.com/brentp/cyvcf2/tarball/276e642245777523b3acd42075a2857da90f1bf3#egg=cyvcf2-0.7.0',
    ],
    install_requires=[
        'cosmos-wfm==2.0.10.indraniel4https',
        'click==6.7',
        'clint==0.5.1',
        'matplotlib==1.5.3',
        'scipy==0.18.1',
        'pandas==0.18.1',
        'seaborn==0.7.1',
        'Cython==0.25.2',
        'cyvcf2==0.7.0',
    ],
    entry_points='''
        [console_scripts]
        yaps2=yaps2.cli:cli
    ''',
    packages=find_packages(exclude=('tests', 'docs')),
    include_package_data=True,
)
