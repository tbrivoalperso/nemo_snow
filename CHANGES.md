*******
Changes
*******

.. todo::

  
# List the main additions of the new 4.2 release - February 2022

This document lists the improvements and new features of the 4.2 release compared to the previous 4.0.
If you are willing to port an existing configuration in order to start using this new release, it is sugggested to also look at the 4.2 Migration Guide (link to add).

Improvements and new features of 4.2 release

## KERNEL

- Quasi Eulerian Coordinates
- New HPG schemes (imrpoveents and new option)
- Preparatory stages for new RK3 time stepping scheme (includes time-pointer changes, DO LOOP macros, changes in main arrays dimensions)
- Full Shallow Water setup
- New vertical scale factors management in time reducing meùory footprint : added key qco (Quasi Eulerian Coordinates), key linssh, see Migration Guide

## AIR SEA INTERACTIONS

- Currents feedbacks
- Mass-flux convection scheme (still not compatible with TOP & PISCES)
- Bulk improvements 
- Atmospheric Boundary layer model (1D vertical as for now) & tools improvements
- Wave forcing improvements

## AGRIF zooms

- Improvement of Agrif for global configurations (periodic, north fold zoom, HPC), 
- Allow AGRIF for multiple vertical grids
- Updated Nesting Tools to set up the AGRIF zooms

## Sea-ice SI3

- EAP & VP rheology (V&V to complete)
- Melt ponds (preliminary implementation)

## ENHANCEMENTS

- Ocean column properties in NEMO-ICB
- OSMOSIS
- Update internal tidal mixing 

## Tracers and biogeochemistry TOP

- Ice sheet iron sources
- Scheme for vertical penetration of visible light

## DATA interface

- OBS modifications

## Input Output Manager

- Use XIOS to read & write restart file (allowing to produce a unique restart file while using domain decomposition)

## High Performance Computing HPC


- MPI Communication cleanup & improvements using MPI3
- Reduce memory footprint
- Improved computational performance of solar penetration scheme
So as some new features implemented, and now needing further work to produce actual performance improvements:
- Extra Halo extension
- Mixed precision preparatory phase
- Loop fusion
- Tiling


## VERIFICATION & VALIDATION

- Tests cases now available with most  new developments
- SETTE validation script improvements
- OASIS test case for ocean atmosphere coupled interface
