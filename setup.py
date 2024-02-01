try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name='mgpgtools',
    version='0.0.1',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=[
        'extern',
    ],
    install_requires=[
        'gfapy',
        'biopython',
        'pandas',
        'toytree',
        'prettytable'],
)
