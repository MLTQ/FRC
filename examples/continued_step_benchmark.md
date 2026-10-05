# Checked-continuation performance protocol

Uses the complete immutable step recorder with ContinuedModalGauss and an
explicitly extended source closure. The temporary protocol binding is restored
on success or failure; use serially in one process. Physics and tolerance
arguments stay identical. All actual stages/accepted states are saved, and
failures are terminal with their last accepted state. This is a measured solve
control, not a full physical trajectory audit or formation result.
