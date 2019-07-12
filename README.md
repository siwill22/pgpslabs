# pgpslabs

Code to make simple reconstructions of slab geometries using gplates-format topological plate reconstructions. 

The idea is principally inspired by various papers of Derek Thorkelson to map the extent of slab windows through geological time. The windows form where mid-ocean ridges intersect with subduction zones. This code maps points in the slabs themselves, so that slab-windows are visualised as gaps between different slab segments. 

Assumptions:
- all calculations assume a single dip angle throughout all slabs
- convergence rates are taken from the plate reconstruction

In current version, the slabs are represented as points along lines of equal subduction age (currently no method is implemented to turn these into a surface). The attributes at each point are the age of subuction, and optionally the age of the seafloor (requires seafloor age grids from which values would be interpolated).
