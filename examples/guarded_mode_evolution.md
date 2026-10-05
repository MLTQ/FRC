# Complete trajectory with bounded proposal work

Wraps the frozen accepted-step trajectory protocol with GuardedModalGauss,
adding its complete source closure explicitly. All original input validation,
stage and endpoint checks, immutable accepted-step arrays/hashes, progress,
terminal failure and final states remain unchanged. Uses temporary serial
bindings restored on any exit. A separate process is required for parallel
variants. No physics or acceptance tolerance is changed by this entry point.
