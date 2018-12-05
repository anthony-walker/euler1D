*********************
Shock Tube Calculator
*********************

Description
===========

This module calculates analytical solutions for a Sod shock tube problem. 

A more detailed description can be found here_.

.. _here: https://gitlab.com/fantaz/simple_shock_tube_calculator 


Usage
=====

.. code:: python

  from shocktubecalc import sod

  # left_state and right_state set p, rho and u
  # geometry sets left boundary on 0., right boundary on 1
  # and initial position of the shock xi on 0.5
  # t is the time evolution for which positions and states in tube should be calculated
  # gamma denotes specific heat
  # npts is number of points to be calculated for rarefaction wave
  # Note that gamma and npts are default parameters (1.4 and 500) in solve function.

  positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.),
                                         geometry=(0., 1., 0.5), t=0.2, gamma=1.4, npts=500)








