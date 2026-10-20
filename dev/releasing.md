# Releasing

The package version comes from git tags ([hatch-vcs](https://github.com/ofek/hatch-vcs)),
so there is no version to bump in `pyproject.toml`. On a tag `v0.2` the version is
`0.2`; between tags it is a dev version such as `0.3.dev4+g1a2b3c4`.

To cut a release, tag `main` and push the tag:

```bash
git switch main && git pull
git tag v0.2
git push origin v0.2
```

The [release workflow](../.github/workflows/release.yml) checks that the package
version matches the tag, runs the tests and publishes a GitHub Release with
generated notes. Tags are `v<major>.<minor>[.<patch>]`.

## Caveats

- **Stale version after a new tag.** pixi installs gaia-cli in editable mode, and it
  doesn't rebuild that install when only the tag changes (it has no equivalent of uv's
  `cache-keys`). Run `pixi reinstall gaia-cli` to refresh what
  `importlib.metadata.version("gaia-cli")` reports. Provenance is unaffected in a
  checkout: `gaia:cli-version` records the git commit there (see
  [provenance-metadata.md](provenance-metadata.md)).
- **CI needs the tags.** `actions/checkout` fetches a single commit by default, which
  gives a meaningless version. The workflows set `fetch-depth: 0`.
- **The lock file doesn't record the version.** `pixi.lock` (format v7) lists the local
  package without a version, so new commits and tags don't change it.
