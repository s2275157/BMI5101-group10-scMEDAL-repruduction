# Results policy

Commit only small, reviewable results such as aggregated CSV/JSON tables and selected figures.

Do not commit:

- AML expression matrices or annotations redistributed from the authors.
- `.npy`, `.h5ad`, checkpoints, or trained weights.
- PBS output logs containing unnecessary environment or host information.
- credentials, SSH keys, tokens, or private paths beyond the documented reproducibility manifest.

Canonical large outputs remain on NSCC scratch and must be backed up separately according to NSCC retention policy.

