# Current/thermal protocol driver

Runs either the full explicit protocol, a kind (original/wave/symbol) or one
exact named case. Independent groups write disjoint files and do not claim a
complete collection. An unknown group/case is rejected explicitly.
Assembly requires every declared configuration, matching
frozen source, NPZ hash and explicit successful/rejected status, plus the full
principal-symbol record. Missing or rejected cases cannot masquerade as success.

Source files must remain frozen while any measurement runs. A protocol-finished
collection can include scientifically informative rejections; each retains its
accepted time and actual failure. This is a bounded verification protocol, not
an automatic formation or physical-model acceptance gate.
