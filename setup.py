from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent


setup(
    name="mlva-seer",
    version="0.1.0",
    description="MLVA/VNTR genotyping from sequencing reads and assemblies",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url="https://github.com/microbemarsh/mlva_seer",
    project_urls={
        "Repository": "https://github.com/microbemarsh/mlva_seer",
        "Issues": "https://github.com/microbemarsh/mlva_seer/issues",
    },
    author="MLVA Seer contributors",
    license="GPL-3.0-only",
    keywords=["mlva", "vntr", "amplicon", "genotyping", "bioinformatics"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
    packages=find_packages(include=["mlva_seer", "mlva_seer.*"]),
    install_requires=[
        "numpy>=1.24",
        "parasail>=1.3.4",
        "pysam>=0.22",
        "sassy-rs>=0.2.4",
    ],
    extras_require={
        "dev": ["pytest"],
    },
    entry_points={"console_scripts": ["mlva-seer=mlva_seer.cli:main"]},
)
