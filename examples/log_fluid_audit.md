# Independent saved-state and bitwise replay audit

Requires the exact frozen source identity and all raw-array hashes. Recreates
initial state/rate and each complete case's last RK4 step bit for bit. A separate
order18 root-product/Legendre sampler measures all saved snapshots, local true
physical energy and initial physical energy rate. Saved history totals replay
at every snapshot. Reports rates relative to their signal and absolute RMS;
equilibrium residuals must not be disguised by dividing by a large background.
Failed trajectories retain their status and are never counted as completed runs.
