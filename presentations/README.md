# Editable course presentation

## Version 6 — 12 October 2026

[Download scMEDAL_AML_reproduction_bilingual_v6.pptx](scMEDAL_AML_reproduction_bilingual_v6.pptx)

This revision follows the version-5 design and evidence. Every slide has the original Chinese notes followed by an English speaking script. The cover gives the full name Acute Myeloid Leukemia (AML), with a brief explanation in its notes. Small source captions now appear in the notes, while inline reference numbers and the reference slide remain visible. The submission deliverable is the editable PPT. No group-contribution page is included.

The [previous version 5](scMEDAL_AML_reproduction_results_v5.pptx) remains available for revision history.

## Slide map

| Slides | Content |
|---|---|
| 1 | Scientific title and reproduction scope |
| 2–4 | Background, representative integration methods, and donor/disease confounding |
| 5–7 | Mixed-effects motivation, separate FE/RE subnetworks, and objectives |
| 8–9 | Data, compute, and evaluation design |
| 10–11 | Five-fold test ASW figure with 95% t intervals, then mean ± SD table |
| 12–14 | FE/RE UMAPs colored by batch, cell type, and patient group |
| 15 | Donor-swap reconstruction as an as-if simulation |
| 16–18 | Selected gene shifts, five-fold direction stability, SAMSN1 target-donor sensitivity |
| 19–20 | Results discussion, limitations, and prospective experiments |
| 21–23 | Conclusion, GitHub archive, and references |
| 24 | Appendix: configuration and environment |
| 25 | Appendix: disease-analysis literature extension; hidden in slide-show mode |

There are 23 main slides and two appendices. Presenter notes provide interpretation and evidence context. Tables, text, and diagrams remain editable; scientific figure panels are embedded images, with standalone figures and their plotting script in the repository.

## Editing guidance

- Fill in any presenter details and adapt pacing to the required presentation duration.
- For shorter talks, keep slide 11 as backup after explaining slide 10; slides 13/14 can be summarized together.
- Keep the distinction between FE's batch-invariant objective and RE's batch-specific objective.
- Retain scANVI's label-use qualification, the descriptive reconstruction scope, and the repeated-target-donor limitation.
- Weight tuning, donor-held-out testing, combined FE+RE prediction, and encoder extensions are prospective work.

The file was exported, structurally checked, and all 25 slides were rendered for visual review. Native PowerPoint playback should be checked after any subsequent edits.

For figure provenance and fuller discussion, see [PRESENTATION_RESULTS.md](../docs/PRESENTATION_RESULTS.md).
