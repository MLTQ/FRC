# Preserve one complete signed workflow from a stopped batch

Requires an individually saved case_NN.json with status complete, exact equality
to its parent batch record, verified complete raw arrays and unchanged original
calculation sources. An incomplete/failed case cannot be selected. Copies all
original arrays without recomputation or modification.

The new archive explicitly records selection_only, parent_batch_terminal,
parent index and all parent summary/provenance/source hashes. Its terminal
status describes the selected complete signed workflow only; the parent
batch's incomplete status and every unfinished record remain unchanged.

The source archive includes this extraction and the entire original calculation
closure. Fresh extraction replay can verify every byte/non-runtime metric.
This is evidence preservation, not a new physical solve, unit test, acceptance
threshold or fabricated completion of an interrupted calculation.
