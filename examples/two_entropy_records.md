# Joint two-entropy provenance and checkpoint records

Traverses the complete local study/reference/audit/test import closure. Sources
and companion docs are archived before measurement; sorted SHA256 entries
define the source identity. Python/NumPy/SciPy versions and exact raw hashes are
retained. Checkpoint publication uses the unchanged atomic serializer pattern
from face_radial_records: distinct running-step arrays precede metadata, and
terminal records are separate. No older source closure is modified.
