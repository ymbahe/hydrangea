"""Demonstration script to plot the location of subhaloes around a cluster."""

import numpy as np
import hydrangea as hy

sim = hy.objects.Simulation(index=0)

subfind_file = sim.get_subfind_file(29)
snapshot_file = sim.get_snap_file(29)

fof_cl = hy.SplitFile(subfind_file, 'FOF', read_index=0)
subhaloes = hy.SplitFile(subfind_file, 'Subhalo')



