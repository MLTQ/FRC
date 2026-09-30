# Common radial-face provenance and checkpoints

Traverses the actual study/audit/test import closure, preserving all older
source identities. Archives the Python sources and companion documentation.
Records package versions, exact source hashes and raw-array hashes.

Running checkpoints use distinct accepted-step array filenames and publish
their metadata atomically only after the array file is complete. The previous
checkpoint remains recoverable after interruption. Terminal records have a
separate name and never masquerade as a complete running checkpoint. There
is one worker per case; no concurrent writers to the same case are supported.
