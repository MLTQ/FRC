# Complete opt-in tensor-mass physical workflow

Runs the unchanged full direct-metric signed protocol with TensorMassGeometry
substituted only during fixture initialization. Every physical equation, actual
density/mass, boundary, original residual and tolerance remains inherited.
Temporary serial bindings restore on success, failure or interruption.

Explicit source closure includes both geometries, full dispatch and tensor
inverse implementation. Numerical inverse results may differ in rounding;
complete actual inputs and geometry operators must remain identical, and
independent full physical fields/rates and original acceptance checks must be
compared. Each adopted workflow needs a fresh complete bitwise replay. No unit
tests are added or run. Resource feasibility and physical accuracy are separate
gates, and coarse pilot success does not authorize fine-grid adoption.
