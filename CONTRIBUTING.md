# Contributing to skinflame

Thanks for your interest! Contributions of all kinds are welcome:
bug reports, feature requests, documentation fixes, new estimators, and
real-data demos.

## Reporting bugs

Open an issue at
https://github.com/vp-thehuman/skinflame/issues/new and
include:

- The output of `python -c "import skinflame, sys; print(skinflame.__version__, sys.version)"`
- A minimal reproducible example (ideally < 30 lines, using
  `simulate_ad_cohort` so anyone can rerun it).
- The full traceback.

## Asking for help

Two routes:

1. **GitHub Discussions** at the repo, for usage questions and "is
   this the right method for my data?" conversations.
2. **GitHub Issues**, if you think you've found a bug.

Please don't email maintainers privately for usage questions — keeping
discussion public means the next person with the same question finds the
answer.

## Proposing a new feature

Open an issue first describing the feature, why it belongs in
`skinflame` (vs. a downstream notebook), and a sketch of the API. We'll
discuss before you sink time into a PR.

## Development setup

```bash
git clone https://github.com/vp-thehuman/skinflame.git
cd skinflame
pip install -e ".[demo,dev]"
pytest                  # 15 tests should pass
ruff check src tests    # zero warnings expected
```

## Code style

- Python 3.10+, type hints encouraged on public APIs.
- `ruff` for linting and import ordering. Config in `pyproject.toml`.
- Keep the public API in `skinflame/__init__.py` minimal and stable.
- Docstrings: NumPy style. Public functions need a one-line summary,
  parameter table, and a short example or note where it's non-obvious.

## Tests

- New estimators need a unit test that recovers a known ground truth
  on `simulate_ad_cohort` data within a stated tolerance.
- New utilities need a unit test for the obvious edge cases.
- Reviewers run `pytest` before merging. If your change makes the suite
  slower than ~5 seconds total, mark slow tests with
  `@pytest.mark.slow` and add a `pytest.ini` marker.

## Pull requests

- Fork → branch → PR.
- Reference the issue your PR closes.
- Keep PRs focused. One feature or one fix per PR.
- CI must be green before merge.

## Code of conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
By participating you agree to abide by its terms.

## License

By contributing you agree that your contributions will be licensed
under the MIT License (see [LICENSE](LICENSE)).
