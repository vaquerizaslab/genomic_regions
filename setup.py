import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('genomic_regions/version.py').read())


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


setup(
    name='genomic_regions',
    version=__version__,
    description='Consistently handle genomic regions',
    packages=find_packages(exclude=["test"]),
    install_requires=[
        'pandas',
        'numpy',
        'pybedtools>=0.8.0',
        'pysam',
        'future',
        'pyBigWig',
        'intervaltree',
    ],
    scripts=['bin/convert-regions'],
    cmdclass={
        'clean': CleanCommand
    }
)
