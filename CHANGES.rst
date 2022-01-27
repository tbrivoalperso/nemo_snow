*******
Changes
*******

.. todo::

   List the main additions of the new 4.2 release - February 2022

   This document lsts the improvements and new features of the 4.2 release compared to the previous 4.0.
   If you are willing to port an existing configuration in order to start using this new release, it is sugggested to also look at the 4.2 Migration Guide (link to add).

Improvements and new features of 4.2 release

## KERNEL

- Quasi Eulerian Coordinates

- New vertical scale factors management

- New HPG schemes, Runge Kutta (preparatory stages includes time-pointer changes, DO LOOP macros

- Full Shallow Water equations setup

## AIR SEA INTERACTIONS

- Currents feedbacks

- Mass-flux convection scheme (still not compatible with TOP & PISCES)

- Bulk improvements 

- ABL1D model & tools improvements

- Wave forcing improvements

- Atmospheric boundary layer

## AGRIF

- Improvement of Agrif for global configurations (periodic, north fold zoom, HPC), AGRIF vertical interpolation

## SI3

- EAP & VP rheology (V&V to complete)

- Melt ponds (preliminary implementation)

## ENHANCE

- Ocean column properties in NEMO-ICB

- OSMOSIS

- Update tidal mixing 

## TOP

- Ice sheet iron sources

- Scheme for vertical penetration of visible light

## DATA interface

- OBS modifications

## IO Manager

- Restart read write with XIOS

## VALIDATION

- Tests cases now available with most  new developments

- SETTE validation script improvements

- OASIS test case for ocean atmosphere coupled interface

## HPC

- Extra Halo extension

- Mixed precision preparatory phase

- MPI Communication cleanup & improvements using MPI3

- Loop fusion

- Tiling

- Improved computational performance of solar penetration scheme
