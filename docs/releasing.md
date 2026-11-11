# Release checks

Use this guide after the mathematical tests and user documentation agree with the release source. Keep the source revision, dependency versions, test results and distribution hashes together so another person can reproduce the checks.

## Check the release inputs

For alggeom 0.1.0, the core dependency floor is algroots 0.4.0; the optional real profile requires semialg 1.5.0. Release those backend versions before testing an index-only installation of the corresponding profile. Source checkout tests are useful during development, but do not establish that users can resolve the published requirements.

Check that [package metadata](../pyproject.toml), [backend inputs](../ci/backends.toml), the [README](../README.md) and [release notes](../RELEASE_NOTES.md) describe the same capabilities and requirements. The metadata declares `GPL-3.0-only` and includes the root [LICENSE](../LICENSE).

## Test the dependency profiles

Use separate clean environments for core and real tests. The core environment must not contain semialg:

```sh
python tools/install_backends.py --published --profile core
python -m pip install ".[test]"
python tools/check_backends.py --profile core --ground-types flint
python -m pytest --without-real
```

In the real environment:

```sh
python tools/install_backends.py --published --profile real
python -m pip install ".[test,real]"
python tools/check_backends.py --profile real --ground-types flint
python -m pytest
```

Run the [CI matrix](../.github/workflows/tests.yml) across its configured operating systems, Python interpreters, minimum/latest arithmetic dependencies and ground coefficient backends. Inspect the results for the release revision. The [testing guide](testing.md) explains the independent expectations, optional oracles and skip reasons.

## Validate the documentation

```sh
ruff check .
ruff format --check .
python tools/generate_api.py --check
```

The suite executes the README, Python tutorial blocks and catalog examples. It also resolves documentation links into the source tree. Check mathematical results and stated limits as well as successful execution. Generated API documentation should match the public source.

## Build and install the distributions

In the real environment, install the build tools and create both formats:

```sh
python -m pip install build twine
python -m build
python -m twine check dist/*
python -m pip install --force-reinstall --no-deps --no-index --find-links dist alggeom
python tools/check_installation.py --profile real
```

The installation checker rejects imports from the package or sibling backend source checkouts, runs each advertised application in a fresh process and executes the documented Python blocks. It retains assertions and per-workflow timeouts.

Rebuild a wheel from the source distribution, then check that installed wheel too:

```sh
python -m pip wheel --no-deps dist/alggeom-0.1.0.tar.gz --wheel-dir rebuilt
python -m pip install --force-reinstall --no-deps --no-index --find-links rebuilt alggeom
python tools/check_installation.py --profile real
```

Confirm that the source distribution includes documentation, examples, tests, backend requirements and build tooling. Check that the wheel includes the license and that its metadata preserves the optional semialg requirement. Repeat installed core checks in the semialg-free environment.

## Record the outcome

Keep full test counts, failures, skips, dependency versions, installed package paths and artifact hashes. Record unavailable external oracles and unrun platforms separately from completed checks. Performance reports should include input dimensions and timing distributions, with no timing thresholds substituted for mathematical correctness.

## Publish the checked revision

The [publishing workflow](../.github/workflows/publish.yml) builds and tests the candidate before uploading its artifacts. A manual run defaults to `publish=false` for a build without an upload. Publication requires a tag matching the version in `pyproject.toml`; for this release, use `v0.1.0`.

Configure the existing alggeom project's [PyPI Trusted Publisher](https://docs.pypi.org/trusted-publishers/adding-a-publisher/) with owner `BhuvaneshBhatt`, repository `alggeom`, workflow filename `publish.yml` and environment `pypi`. The filename and environment must match the publishing job.

After the checks pass and the tag contains the reviewed source, start publication explicitly:

```sh
gh workflow run publish.yml --repo BhuvaneshBhatt/alggeom --ref v0.1.0 -f publish=true
```

Publishing a GitHub Release also triggers this workflow and its upload. Follow the run through artifact checks, PyPI publication and the indexed installation check. If authentication fails, correct the publisher configuration and rerun the failed publishing job to reuse its successful build artifacts.
