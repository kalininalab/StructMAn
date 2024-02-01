from setuptools import setup, find_packages

path_to_readme = "README.md"

if path_to_readme is not None:
    with open(path_to_readme, "r") as desc_file:
        long_description = desc_file.read()
else:
    long_description = ''

path_to_version_file = "./structman/_version.py"

with open(path_to_version_file) as version_file:
    exec(version_file.read().strip())

setup(
    name="StructMAn",
    version=__version__,
    description="Structural Mutation Annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='LGPL-2.1',
    author="Alexander Gress",
    maintainer="Alexander Gress",

    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    setup_requires=['setuptools_scm'],
    include_package_data=True,
    install_requires=[
        "biopython>=1.79",
        "python-igraph==0.10.8",
        "matplotlib>=3.7.1",
        "numpy>=1.22.3",
        "psutil>=5.8.0",
        "pymysql>=1.0.2",
        "ray==2.9.1",
        "msgpack>=1.0.3",
        "zstd>=1.5.2.5",
        "pandas>=1.3.3",
        "autopep8>=1.5.7",
        "scipy==1.11.3",
        "requests==2.28.1",
        "more-itertools>=9.1.0",
        "pycairo>=1.23.0",
        "powerlaw>=1.5",
        "biotite>=0.38.0",
        "markdown>=2.6.9"
    ],

    package_data = {
        "": [
            'structman/lib/rinerator/reduce',
            'structman/lib/rinerator/probe',
            'structman/lib/rinerator/reduce_wwPDB_het_dict.txt',
            'structman/resources/pdbba.fasta.gz'
            ]
    },
    python_requires=">=3.8, <4",
    keywords="bioinformatics",
    entry_points={
        "console_scripts": ["structman = structman.structman_main:structman_cli"],
    },
)
