try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup


setup(
    name="mgpgtools",
    version="0.0.1",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "gfapy",
        "biopython",
        "pandas",
        "toytree",
        "prettytable",
    ],
    entry_points={"console_scripts": ["mgpgtools = main:main"]},
)
