# Explicit disk-failure case recovery

Recomputes only explicitly requested protocol indices whose case files are
absent, under the unchanged frozen scientific source. It refuses existing case
files or an existing recovery record and records the storage-related reason.
This is external job recovery, not a hidden numerical time-step retry or a
pressure/density repair. Completed original cases remain untouched.

An explicit assembly operation verifies all ten case identities/raw hashes and
creates the missing complete protocol only when every case is complete. The
original worker must be terminal before assembly; this helper does not stop or
replace live work. Its source lies outside the frozen calculation closure and
is recorded separately in the recovery/checkpoint evidence.
