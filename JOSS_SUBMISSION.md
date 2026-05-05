# JOSS submission — step-by-step

This file walks you through getting `skinflame` accepted at the **Journal
of Open Source Software** (JOSS). Tick the checkboxes as you go.

JOSS is unusual among journals: the entire review is on a public GitHub
issue, the editor is a volunteer scientist, the reviewers are real
researchers who *use* your software, and the bot (`@editorialbot`) drives
the workflow. Median time from submission to publication is 90 days.

---

## Phase 0 — Decide if `skinflame` qualifies (5 min)

JOSS's gatekeeping rules. All must be true:

- [x] **Open source**, OSI-approved license. (MIT — done.)
- [x] Software is **research software**, not just a generic library. AD
      multi-omic mediation is a research topic; this passes.
- [x] **Substantial scholarly effort**. JOSS rejects roughly a quarter of
      submissions at the *pre-review* stage for being too small. The bar
      is informally "more than a weekend project". `skinflame` ships
      ~2,500 LOC, 3 distinct estimators (classic / HIMA / pathway),
      DAG-aware adjustment, a Streamlit demo on real GEO data, and 15
      unit tests. This is comfortably above the bar.
- [x] **Author commitment**. JOSS expects you'll respond to issues for
      ~5 years. Realistic? If yes, proceed.
- [ ] **You are not on the JOSS editorial board.** (You're not.)

If any of the above is iffy, file a pre-submission inquiry first:
https://github.com/openjournals/joss/issues/new?template=pre-review-question.md

---

## Phase 1 — Make the repo JOSS-ready (one evening)

### 1.1 Personal info — search & replace these strings

In `pyproject.toml`, `LICENSE`, `CITATION.cff`, `paper.md`, `README.md`:

| Placeholder                                | Replace with                          |
| ------------------------------------------ | ------------------------------------- |
| `REPLACE_WITH_YOUR_NAME`                   | Your full name                        |
| `REPLACE_WITH_YOUR_LAST_NAME`              | Surname                               |
| `REPLACE_WITH_YOUR_FIRST_NAME`             | Given name                            |
| `REPLACE_WITH_YOUR_AFFILIATION`            | University / institute / "Independent researcher" |
| `REPLACE_WITH_YOUR_GITHUB`                 | GitHub username                       |
| `0000-0000-0000-0000` (in CITATION.cff and paper.md) | Your real ORCID iD          |
| `your-org/skinflame` (URLs in pyproject)   | `<your-github>/skinflame`             |

If you don't have an ORCID iD: register one in 2 minutes at
https://orcid.org. JOSS *requires* it for the corresponding author.

### 1.2 Repository hygiene JOSS reviewers always check

Reviewers walk through a fixed checklist (you can preview it here:
https://joss.readthedocs.io/en/latest/review_checklist.html). The
non-obvious items:

- [x] **Statement of Need** is in `paper.md` (already there).
- [ ] **Community guidelines** in repo root or `docs/`:
      - [x] `CONTRIBUTING.md` (added by this scaffold — see file)
      - [x] `CODE_OF_CONDUCT.md` (added — Contributor Covenant 2.1)
      - [ ] How to file a bug → covered in `CONTRIBUTING.md`
      - [ ] How to seek support → covered in `CONTRIBUTING.md`
- [x] **Installation instructions**: `pip install skinflame` (README + `docs/quickstart.md`).
- [x] **Example usage**: README has a runnable example; `docs/quickstart.md`
      has more.
- [x] **Functionality documentation**: docstrings on all public functions
      + `docs/theory.md`.
- [x] **Tests**: `pytest tests/` (15/15 passing).
- [x] **Automated tests run in CI**: `.github/workflows/ci.yml`.
- [x] **API documentation**: docstrings (you can optionally add MkDocs
      or Sphinx, but it's not required).

### 1.3 Push to GitHub

```bash
cd skinflame
git init
git add -A
git commit -m "Initial release: skinflame v0.1.0"
git branch -M main
git remote add origin git@github.com:<your-github>/skinflame.git
git push -u origin main
```

Make sure GitHub Actions runs and stays green on `main`. The CI badge
on your README is what reviewers click first.

### 1.4 Tag a release

```bash
git tag -a v0.1.0 -m "skinflame v0.1.0 — JOSS submission"
git push --tags
```

Then on GitHub: **Releases → Draft a new release → choose tag v0.1.0 →
Publish release.**

### 1.5 Get a Zenodo DOI for the release

1. Sign in to https://zenodo.org with your GitHub account.
2. **Settings → GitHub → toggle `skinflame` to On.**
3. Go back to GitHub and *re-publish* the v0.1.0 release (or push a
   v0.1.1 tag). Zenodo will archive it and mint a DOI within ~5 minutes.
4. Copy the DOI badge MD snippet from the Zenodo "DOI" tab and paste it
   at the top of your README.

You now have a citable DOI like `10.5281/zenodo.12345678`.

### 1.6 Final sanity check

```bash
pip install -e ".[demo,dev]"
pytest                              # all green
ruff check src tests                # zero warnings
streamlit run demo/streamlit_app.py # demo works in your browser
```

If any of these fails, fix it now — JOSS reviewers will run the same
commands.

---

## Phase 2 — Write the paper (2-4 hours)

`paper.md` is already drafted. JOSS papers are tightly constrained:

- **Length**: 250-1000 words (renders to ~2 typeset pages with figures).
  More is allowed but reviewers will push back.
- **Required sections** (you have all of them):
  1. *Summary* (high-level, audience = a non-specialist scientist)
  2. *Statement of Need* (why this exists, who needs it, what's missing
     from existing tools — this is the section reviewers scrutinise most)
  3. (Optional but expected) *Methods* / *Functionality*
  4. *References* — must use BibTeX, must include DOIs where they exist
- **Author block**: each author needs a name, ORCID, and affiliation.
- **Figures**: optional. If you add any, put PNG/PDF in `paper/` and
  reference with `![caption](paper/figure.png)`.

### Tighten paper.md before submission

- [ ] Read it aloud. Cut anything that doesn't help a reviewer answer
      "what does this software do, who is it for, why is it needed?"
- [ ] Confirm every cited paper has a DOI in `paper.bib`.
- [ ] Add 1-2 sentences in *Statement of Need* explicitly naming the
      Python tools `skinflame` competes with (or notes the absence of
      Python competitors). E.g.: *"While R offers `mediation`, `HIMA`,
      and `bama`, no Python equivalent bundles classical, high-dimensional,
      and pathway-level mediation behind a single API."*
- [ ] Run a local render to check it compiles. JOSS uses
      [Open Journals' `inara` Docker image](https://github.com/openjournals/inara):

      ```bash
      docker run --rm -it \
        -v $PWD:/data \
        -u $(id -u):$(id -g) \
        openjournals/inara
      ```

      Output PDF lands at `paper.pdf`. Reviewers run this same command.

---

## Phase 3 — Submit (15 min)

1. Go to https://joss.theoj.org/papers/new
2. Sign in with GitHub.
3. Fill in the form:
   - **Repository URL**: `https://github.com/<your-github>/skinflame`
   - **Software version**: `v0.1.0`
   - **Branch**: `main` (or whichever holds `paper.md`)
   - **Submission Track**: pick `Computer science` *or*
     `Health, Life & Earth Sciences` — the latter usually finds you
     reviewers who know AD/genomics. If unsure, pick `Health, Life & Earth Sciences`.
   - **Author / Co-author**: just you, unless others contributed.
4. Click **Submit**. The bot opens an issue at
   `openjournals/joss-reviews` titled "[REVIEW]: skinflame: ..."
   within a minute.

---

## Phase 4 — Pre-review (1-2 weeks)

The JOSS Editor-in-Chief assigns a **handling editor**. The handling
editor and `@editorialbot` will run pre-review checks:

- License OK? (MIT ✓)
- Substantial scholarly effort? (You may need to argue this — see below)
- Statement of Need present? (✓)
- Tests present? (✓)
- Authors all have ORCIDs? (Make sure!)
- Repo public and active? (✓)

If pre-review passes, the editor recruits 2 reviewers (it can take 2-4
weeks to find them). If pre-review fails, the editor explains why in the
issue and gives you a chance to fix it.

**If asked "is this substantial enough?":** answer with a short comment
listing: LOC count, number of distinct algorithms, number of unit tests,
the GSE121212 real-data demo, and the explicit fact that no Python
equivalent exists for multi-omic mediation. This usually settles it.

---

## Phase 5 — Review (4-8 weeks)

Each reviewer posts a checkbox checklist as a single comment on the
review issue. They tick each item as they verify it. You'll get one of:

- **Box ticked** — no action.
- **Box unticked + comment** — they hit a problem. Your job is to fix it
  (in your repo) and reply on the issue with the commit/PR link.

Common reviewer requests for a package like this:

| Request                                   | Fix                                       |
| ----------------------------------------- | ----------------------------------------- |
| "I get an `ImportError` on Python 3.10"  | Add a CI matrix row for 3.10 (already there) |
| "Tests pass but coverage is low"          | Add 2-3 tests for edge cases              |
| "How do I cite this?"                     | `CITATION.cff` (already added)            |
| "How do I contribute?"                    | `CONTRIBUTING.md` (already added)         |
| "What's the difference vs. R `HIMA`?"     | Add a "Comparison to existing tools" subsection in README or paper.md |
| "Streamlit demo errors when I click X"    | Fix and push                              |
| "Statement of Need is too vague"          | Sharpen the paper sentence-by-sentence    |
| "Add a CHANGELOG.md"                      | Trivial; do it                            |
| "Bump dependency floor on numpy"          | Update `pyproject.toml`                   |

Reply within ~7 days of each comment. Slow responses are the #1 reason
JOSS submissions stall.

---

## Phase 6 — Acceptance & publication (1-2 weeks)

When both reviewers tick all boxes:

1. The editor recommends acceptance to the EiC.
2. `@editorialbot` asks you to do a final release with the *exact*
   version they're going to publish. Tag it (e.g. `v0.1.1`) and push.
3. Zenodo auto-archives the new release; copy the DOI into the review
   issue.
4. The bot generates the typeset PDF. Read it, OK it.
5. The EiC accepts; the bot publishes. Within ~24h your paper is live at
   `https://joss.theoj.org/papers/10.21105/joss.0XXXX` with a real DOI.

You now have a peer-reviewed first-author publication.

---

## CV bullets unlocked at acceptance

```
PUBLICATIONS
  <Last, F.> (2026). skinflame: Multi-omic causal mediation analysis
  for atopic dermatitis and inflammation research. Journal of Open
  Source Software, 11(XX), XXXX. https://doi.org/10.21105/joss.0XXXX

SOFTWARE
  skinflame (2026) — Multi-omic causal mediation library (Python).
  PyPI · GitHub · Zenodo (DOI: 10.5281/zenodo.XXXXXXX) · CI on Linux,
  3.10/3.11/3.12 · MIT.
```

---

## What to do if it gets rejected

JOSS rejection rates are ~10-15%, mostly at pre-review for "scope".
If it happens:

- Read the editor's stated reason carefully. It is almost always one of:
  scope (too small), license (not OSI), or "this is an application not a
  library".
- For `skinflame` the most likely failure mode is the editor judging
  the project too small. Counter-argument: you have multiple distinct
  estimators, a real-data demo, and a research-novel application area.
- If denied, send the same paper to **Bioinformatics Advances** as a
  Software Note (4 pages, ~6 week turnaround) or **BMC Bioinformatics**
  Software article (~10 pages, longer turnaround). Both are still
  peer-reviewed and DOI-bearing.

---

## Realistic timeline

| Week | Action                                    |
| ---- | ----------------------------------------- |
| 1    | Replace placeholders, push to GitHub, tag v0.1.0, get Zenodo DOI |
| 2    | Polish `paper.md`, render with `inara`, submit on JOSS website |
| 3-4  | Editor pre-review; respond to any pre-review comments |
| 5-10 | Review iteration with two reviewers       |
| 11-12| Final release, typesetting, publication   |

That's a published, peer-reviewed paper from "we have a working repo" in
about three months of evenings.
