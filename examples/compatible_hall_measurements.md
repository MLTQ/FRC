# Complete compatible Hall study runner

Runs every explicit spatial or coupled-evolution configuration sequentially.
Each completed case writes its own auditable record and lossless archive; the
collected completed marker is written only after all cases succeed with one
fixed source fingerprint. Progress messages are not a completion marker.
The report must audit coverage/configuration/source/archive identities before
claiming a completed study. This run does not certify FRC formation or the
unimplemented thermal/vacuum/interface physics.
