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

# Stamp the install with the branch and commit it was built from, so every run can name the code
# it is running (see structman/_build_stamp.py for why). Loaded by path, like _version.py: the
# package's dependencies do not exist yet. Never allowed to fail the install.
try:
    import os as _os
    import runpy as _runpy
    _here = _os.path.dirname(_os.path.abspath(__file__))
    _stamp = _runpy.run_path(_os.path.join(_here, 'structman', '_build_stamp.py'))
    _stamp['write'](_os.path.join(_here, 'structman', '_build_info.py'), _stamp['collect'](_here))
except Exception as _stamp_error:
    print(f'Could not stamp the build with its commit: {_stamp_error}')

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
    include_package_data=True,
    install_requires=[
        "biopython>=1.79",
        "python-igraph==0.10.8",
        "matplotlib>=3.7.1",
        # Upper bound is load-bearing: numba 0.66 moved NPDatetime/NPTimedelta
        # out of numba.core.types.scalars, which every numba-cuda release to
        # date (<=0.30.4) still redirects numba.cuda.types.scalars onto. Going
        # to >=0.66 makes `import numba.cuda` -- and therefore rmm, cupy's RMM
        # allocator and all StructGuy GPU training -- fail at import time.
        "numba>=0.60,<0.66",
        "numpy>=1.22.3",
        "psutil>=5.8.0",
        "pymysql>=1.0.2",
        "ray>=2.59.0",
        "msgpack>=1.0.3",
        # base_utils imports both of these -- they are different projects, not
        # two names for the same one: zstd provides the module-level
        # compress/decompress used for column blobs, zstandard provides the
        # ZstdCompressor/ZstdDecompressor streaming API. Declaring only zstd
        # left `import zstandard` failing on a clean install.
        "zstd>=1.5.2.5",
        "zstandard>=0.21.0",
        # Held below 3.0 deliberately: pandas 3 turns on copy-on-write and the
        # str dtype by default, which this code has not been audited against.
        # Lift once that audit happens -- nothing here needs pandas 3.
        "pandas>=2.2.2,<3",
        "autopep8>=1.5.7",
        "scipy>=1.14.1",
        "requests>=2.32.3",
        "more-itertools>=9.1.0",
        "powerlaw>=1.5",
        "biotite>=1.0.1",
        "markdown>=2.6.9",
        "pdfkit==1.0.0",
        "pybind11>=2.13.6",
        "wget==3.2",
        "pydantic>=2.12.3"
    ],

    extras_require={
        # mysqlclient is a C wrapper around libmysqlclient and escapes query
        # parameters far faster than pure-Python pymysql, which dominates the
        # cost of StructMan's very large bulk inserts. It is picked up
        # automatically when importable (see base_utils/mysql_driver.py);
        # pymysql stays as the always-available fallback.
        "fast-mysql": ["mysqlclient>=2.1.0"],
    },

    package_data = {
        "": [
            'structman/lib/rinerator/reduce',
            'structman/lib/rinerator/probe',
            'structman/lib/rinerator/reduce_wwPDB_het_dict.txt',
            'structman/resources/pdbba.fasta.gz',
            'build_smssp/build.sh'
            ]
    },
    python_requires=">=3.8, <4",
    keywords="bioinformatics",
    entry_points={
        "console_scripts": ["structman = structman.structman_main:structman_cli"],
    },
)
