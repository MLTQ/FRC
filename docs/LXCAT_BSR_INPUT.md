# Private LXCat BSR argon input

On 2026-10-16 the user explicitly authorized “Accept for local research.”
The LXCat terms were accepted through the official browser workflow, and the
BSR electron–argon archive was downloaded successfully. This supersedes older
“permission pending” notes; it does not authorize redistribution or commercial use.

Private storage, outside the repository:
`/Users/max/Library/Application Support/FRC/personal-atomic-inputs/lxcat-bsr-argon-2026-10-16/`

The original archive and unchanged text/XML exports are preserved alongside
`provenance.json` and a screenshot of the accepted download workflow. The text
contains 33 process blocks: one elastic momentum-transfer, 31 excitation and one
ionization. The header identifies BSR-500. These are acquired inputs only; state
mapping, energy ranges and suitability for individual reactions remain to be
reviewed before use. Do not infer excited-parent or final-state branching data
from the number of excitation channels. The simulation and visualization are
unchanged, and broader research remains paused.

## Attribution and restrictions

Reference: **BSR database, www.lxcat.net, retrieved on October 16, 2026.**
The provider header additionally requests:

- O. Zatsarinny (2006), “BSR: B-spline atomic R-matrix codes,” Computer Physics
  Communications 174, 273.
- O. Zatsarinny and K. Bartschat (2013), “The B-spline R-matrix method for atomic
  processes: application to atomic structure, electron collisions and
  photoionization,” Journal of Physics B 46, 112001.

Preserve applicable process-specific references when using individual channels.
[LXCat terms](https://us.lxcat.net/instructions/how_reference.php) require database,
site, retrieval date and original-source attribution, and prohibit third-party
distribution and commercial use. Raw tables and reconstructed provider rate
arrays must not be bundled in Git, source archives or shared artifacts.
[Database](https://www.lxcat.net/BSR),
[official download workflow](https://us.lxcat.net/data/download_set.php).

## Exact input identities

| Private file | Bytes | SHA-256 |
|---|---:|---|
| download.zip | 68070 | `0d7058707eeb0ac25ee7d537d6e119bb8a42c9d586a689cb2dc6233e5b0c6123` |
| cross_sections.txt | 183265 | `ed851c2f5f9292ff46ff890069a28b3f4c4edba1068ae130c0a3dd6080758f44` |
| cross_sections.xml | 174279 | `c7ecf9762ccedac10113cc77fddbbf7209e8d4ad2c97cfa4e38241a3b6b4ea3b` |

On resumption, inspect the actual target/product identities and cross-section
energy domains; map only defensible channels into the existing atomic source,
retain explicit data gaps, and compare complete trajectories with independent
baselines. Track this through frc-transport.5.3 and frc-transport.5.4. Do not
start new simulation work solely because this access permission was granted.
