# Bioconda release checklist

This directory stages the recipe that will be copied to
`bioconda/bioconda-recipes/recipes/mlvamaps/meta.yaml`. It is not an installable
channel by itself.

## Current blocker

`mlvamaps` imports `sassy` from the PyPI distribution `sassy-rs>=0.2.6` for its
core primer search. As of 2026-08-27, `sassy-rs` is not available from Bioconda
or conda-forge. Bioconda recipes cannot satisfy runtime dependencies with pip,
so package `sassy-rs` first and have that recipe accepted before submitting
`mlvamaps`.

`spoars`, the other native Python extension imported by the calling pipeline,
is already available from Bioconda. The remaining declared Python packages and
external executables are also available from Bioconda/conda-forge.

## Prepare a release

1. Ensure the version in `mlvamaps/_version.py` matches the version in
   `meta.yaml`; `pyproject.toml` reads the Python version file automatically.
2. Run the test and package checks from the repository root:

   ```bash
   pytest -q
   python -m build
   python -m twine check dist/*
   ```

3. Commit the release, create an annotated tag, and push it:

   ```bash
   git tag -a v0.1.0 -m "mlvamaps 0.1.0"
   git push origin main v0.1.0
   ```

4. Download the immutable tag archive and calculate its checksum:

   ```bash
   curl -L -o mlvamaps-0.1.0.tar.gz \
     https://github.com/microbemarsh/mlvamaps/archive/refs/tags/v0.1.0.tar.gz
   shasum -a 256 mlvamaps-0.1.0.tar.gz
   ```

5. Replace `REPLACE_WITH_RELEASE_TARBALL_SHA256` in `meta.yaml` with that hash.
   The placeholder is intentional until the tag exists; a checksum made from a
   mutable branch archive is not release-safe.

## Submit and validate

1. Fork and clone `bioconda/bioconda-recipes`.
2. Create `recipes/mlvamaps/` and copy in `meta.yaml`.
3. Confirm no package named `mlvamaps` already exists in Bioconda or
   conda-forge.
4. Run the current `bioconda-utils` lint/build workflow from the recipes
   checkout, following the
   [Bioconda contributor documentation](https://bioconda.github.io/contributor/).
5. Open a pull request and disclose that `sassy-rs` is a required native
   dependency with its own prerequisite recipe.

The recipe is `noarch: python` because the `mlvamaps` distribution itself is
pure Python. Platform-specific code is supplied through conda dependencies.
Tests intentionally use imports and CLI help/version commands so they validate
installation without downloading genomes or running expensive analyses.