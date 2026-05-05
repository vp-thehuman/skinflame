# Where to publish `skinflame` — a CV-oriented playbook

Goal: get `skinflame` cited by other researchers, indexed on PyPI, archived
with a DOI, and listed on your CV under both **Software** and
**Publications**. The plan below is ordered roughly by effort vs. payoff,
and the milestones build on each other.

## 0. Pre-flight (1 evening)

- [ ] Create the repo on GitHub: `github.com/<you>/skinflame`. Push everything.
- [ ] Update `pyproject.toml` URLs (`Homepage`, `Issues`) to point at the repo.
- [ ] Replace `"skinflame contributors"` in `pyproject.toml` and `LICENSE`
      with your real name + ORCID iD.
- [ ] Enable GitHub Actions (CI workflow already at `.github/workflows/ci.yml`).
- [ ] Add a one-paragraph "Citing this work" block to the README pointing at
      the Zenodo DOI you'll create in step 2.
- [ ] Tag a release: `git tag v0.1.0 && git push --tags`.

## 1. PyPI — *gives you* `pip install skinflame`

Why: a real PyPI listing is the single most credible signal that a research
software project exists. It also makes citation trivial.

```bash
pip install build twine
python -m build           # produces dist/*.whl and dist/*.tar.gz
twine check dist/*
twine upload dist/*       # needs a PyPI account + API token
```

CV bullet you unlock:
> *Author of `skinflame` (PyPI), an open-source Python library for
> multi-omic causal mediation analysis applied to atopic dermatitis and
> inflammation research.*

## 2. Zenodo — *gives you a DOI you can cite*

Why: Zenodo mints a permanent DOI for each GitHub release. JOSS, bioRxiv,
and most journals require a DOI for software citations.

1. Log in to https://zenodo.org with your GitHub account.
2. Toggle the `skinflame` repo "On" under *Settings → GitHub*.
3. On GitHub, *Releases → Draft a new release* on the `v0.1.0` tag. Zenodo
   archives it automatically and gives you a DOI like `10.5281/zenodo.XXXXX`.
4. Add the DOI badge + bibtex snippet to the README.

CV bullet you unlock:
> *DOI: 10.5281/zenodo.XXXXX (Zenodo).*

## 3. JOSS — *the Journal of Open Source Software* (highest ROI for CV)

Why: JOSS is a peer-reviewed, DOAJ-indexed journal specifically for
research software. A JOSS paper is **3–4 pages**, the review is on GitHub
and is constructive, the typical turnaround is 4–10 weeks, and acceptance
gives you a *citable peer-reviewed publication* with minimal extra writing.

Eligibility checklist (`skinflame` already passes most of these):
- [x] Open-source license (MIT — done).
- [x] Documented installation + minimal example (README + `docs/quickstart.md`).
- [x] Automated tests (`pytest` — done).
- [x] CI on every PR (GitHub Actions — done).
- [ ] Statement of Need (1 paragraph) explaining who the software is for and
      why existing tools aren't enough — write this in `paper.md`.
- [ ] `paper.md` (Markdown) + `paper.bib` (BibTeX) at repo root.
- [ ] At least one substantive author and a "for-real" use case demonstrated
      in the docs (the GSE121212 example covers this).

Submit at: https://joss.theoj.org/papers/new

CV bullet you unlock:
> *<Last>, <First>. (2026). skinflame: Multi-omic causal mediation analysis
> for atopic dermatitis and inflammation research. Journal of Open Source
> Software, 11(XX), XXXX. https://doi.org/10.21105/joss.XXXXX*

## 4. bioRxiv preprint — *gives you a citable methods paper*

Why: A 6–10 page methods preprint on bioRxiv reaches the AD / inflammation /
multi-omics audience directly and gets indexed by Google Scholar within
days. It's the standard route for life-science software in 2026.

Suggested title:
> *skinflame: A Python toolkit for multi-omic causal mediation analysis in
> atopic dermatitis and inflammation*

Sections (tight, ~6 pages):

  1. *Background.* Why mediation? Why multi-omic? Why AD specifically?
  2. *Methods.* Reuse `docs/theory.md`. Add a benchmark vs. R `HIMA` and
     R `mediation` on simulated data (matched estimates within 5%).
  3. *Application to GSE121212.* Show the top FDR-significant
     transcriptomic mediators of lesional vs. non-lesional status.
  4. *Software availability.* PyPI + Zenodo DOI + GitHub.
  5. *Author contributions, funding, conflicts.*

Submit to bioRxiv → **Bioinformatics** subject area. Dual-submit to a peer
journal afterwards (see step 5).

CV bullet you unlock:
> *Preprint on bioRxiv, doi:10.1101/2026.MM.DD.XXXXXX*.

## 5. Peer-reviewed venue (pick one)

After the preprint is up, send the same paper (lightly reformatted) to one
of these. They are ordered by audience fit, not prestige:

| Journal                        | Why it fits                                         | Format |
| ------------------------------ | --------------------------------------------------- | ------ |
| **Bioinformatics** (OUP)       | Standard home for methods/software in genomics      | "Software paper", ~4 pages |
| **Bioinformatics Advances**    | Faster, OA sibling of *Bioinformatics*              | Same   |
| **BMC Bioinformatics**         | Also accepts software papers                        | ~10 pages |
| **PLOS Computational Biology** | Better fit if you push the AD↔mood angle hard       | Full paper |
| **GigaScience**                | Strong for multi-omics + reproducibility            | Full paper |
| **Briefings in Bioinformatics**| Higher impact, prefers reviews + methods together   | Full paper |

For a v0.1 release backed by a JOSS paper, *Bioinformatics Advances* or
*BMC Bioinformatics* is the lowest-friction route to "first-author
peer-reviewed publication". If you also write a benchmark + AD↔mood
discussion section, *PLOS Comp Bio* is realistic.

## 6. Conferences / posters (for visibility, not the paper itself)

- **ISMB / ECCB** (Bioinformatics) — software demo track.
- **RECOMB / RECOMB-Seq** — methods venue.
- **Society for Investigative Dermatology (SID)** annual meeting — direct
  AD audience; abstract deadlines are usually in November.
- **EAACI** (European Academy of Allergy & Clinical Immunology) — strong
  AD ↔ inflammation overlap.
- **PyData / SciPy** — Python software talk; great for the CV "Talks"
  section even if not peer-reviewed.

## 7. Citations + indexing nice-to-haves

- **CITATION.cff** in repo root — GitHub auto-renders it as a "Cite this
  repository" button. (One file, 15 lines, big perceived-rigour win.)
- Submit to **scicrunch.org** / **bio.tools** registries — indexes you
  alongside other bioinformatics tools.
- Add to **Awesome-Bioinformatics** and **Awesome-Causal-Inference**
  lists on GitHub via PR.

## 8. CV scaffolding

Once steps 1-3 are done, your CV should have something like:

> **Software**
> - **skinflame** (2026). Multi-omic causal mediation analysis library
>   (Python). PyPI · GitHub · Zenodo DOI · 15 unit tests · CI on Linux.
>
> **Publications**
> - <Last, First>. *skinflame: Multi-omic causal mediation analysis for
>   atopic dermatitis and inflammation.* Journal of Open Source Software,
>   2026. DOI: 10.21105/joss.XXXXX.
> - <Last, First>. *skinflame: A Python toolkit for ...* bioRxiv 2026.
>   DOI: 10.1101/2026.MM.DD.XXXXXX.

## Suggested timeline (8 weeks, evenings/weekends)

| Week  | Action                                                     |
| ----- | ---------------------------------------------------------- |
| 1     | Polish README, add `paper.md`/`paper.bib`, push to GitHub. |
| 2     | Push to PyPI; cut Zenodo release; add CITATION.cff.        |
| 3     | Write the JOSS paper. Submit.                              |
| 4     | Run benchmark vs. R HIMA + R mediation; write methods paper.|
| 5     | bioRxiv preprint live.                                     |
| 6     | Submit to Bioinformatics Advances or BMC Bioinformatics.   |
| 7-8   | Address JOSS reviewer comments; iterate on v0.2.           |

## Realistic CV impact

- **JOSS paper accepted** ≈ 1 first-author peer-reviewed pub, indexed in
  DOAJ, citable. Strong signal for industry research roles, ML/health
  PhD applications, and bioinformatics positions.
- **bioRxiv + a journal paper** ≈ a second first-author publication and
  Google Scholar visibility for the topic.
- **PyPI + Zenodo DOI** ≈ a "Software" CV section with a stable identifier
  reviewers can click.

If you only do *one* thing from this list, do **JOSS**. It is the highest
CV-impact-per-hour among scientific publishing options for a software
project at this scope.
