from setuptools import setup

setup(
    name="locus_processing",
    version="0.0.4",
    description="Tools for working with locus definition files",
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    url="https://github.com/LUMC/locus_processing",
    platforms=['any'],
    packages=["locus_processing"],
    install_requires=[
        "click>=6.7",
        "marshmallow==2.13.5",
        "requests>=2.18.1",
        "pyyaml>=3.12"
    ],
    tests_requires=['pytest'],
    entry_points={
        "console_scripts": [
            "locus2bed = locus_processing.cli:locus_to_bed",
            "validate_locus = locus_processing.cli:validate_locus",
            "complete_locus = locus_processing.cli:complete_locus"
        ]
    },
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        "License :: OSI Approved :: MIT License",
    ],
    keywords='bioinformatics'
)
