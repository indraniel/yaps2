# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('yaps2/version.py') as f:
    exec(f.read())

setup(
    name='yaps',
    version=__version__,
    description='Yet Another Pipeline Setup',
    long_description=readme,
    author='Indraniel Das',
    author_email='idas@wustl.edu',
    license=license,
    url='https://github.com/indraniel/yaps2',
    install_requires=[
        'git+https://github.com/indraniel/COSMOS2.git@enable-lsf',
        'click',
        'clint',
    ],
    entry_points='''
        [console_scripts]
        yaps2=yaps2.cli:cli
    ''',
    packages=find_packages(exclude=('tests', 'docs')),
    package_data={
        '': ['*.md', 'LICENSE'],
        'yaps2' : ['resources/*'],
    },
)
