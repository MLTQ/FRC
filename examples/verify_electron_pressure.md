# Electron-pressure report

Audits all completed or explicitly expected-failure cases before comparing
time, mesh, proxy-scale and nonlinear-tolerance controls. Writes separately
fingerprinted report JSON and two scientific figures: convergence and actual
pressure-battery field/thermal exchange. Figure labels explicitly retain fixed
density, dual nodal pressure, cell-center Btheta and absent bulk/Hall evolution.
Physical continuous pressure errors are independently integrated and plotted
beside nodal errors. Periodic scalar plots include the closing half dual cell.
Two separately recorded report controls increase the strict failed solve's
iteration budget to 160/320 and verify unchanged inputs and zero accepted time.
These controls share the frozen engine but belong to report-source provenance.
Report completion is not formation completion. Visual inspection is required.
