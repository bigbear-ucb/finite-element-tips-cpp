This folder contains the modification to deal.II source files for purpose of acceleration on assembly.

The original assembly requres:
N_cell x N_local_dof x N_local_dof x N_q multiplications

If we pre-assembly the multiplications in the reference cell, we only have:
N_cell x N_local_dof x N_q multiplications

for 2D bilinear element, the price is scaled down by 4; for 2D biquadratic element, the price is scaled down by 9

The effects of acceleration is specially significant in 3D calculations.
