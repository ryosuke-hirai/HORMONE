# HORMONE
![Test suite](https://github.com/ryosuke-hirai/HORMONE/actions/workflows/test.yml/badge.svg)

Repository for hydrodynamics code HORMONE

## About
HORMONE is a 3D grid-based Godunov-type magnetohydrodynamics code that uses the HLLD approximate Riemann solver.
The main feature is the "hyperbolic self-gravity" solver, which enhances the computational cost for Newtonian self-gravity simulations by a significant margin.
Currently, it can handle Cartesian/cylindrical/spherical coordinates.

## Documentation
Under construction...

## Reference

Hyperbolic self-gravity and original code <br>
![Hirai, Nagakura, Okawa, Fujisawa, 2016, PhRvD, 93, 083006](https://ui.adsabs.harvard.edu/abs/2016PhRvD..93h3006H/abstract)

Equation of state with recombination <br>
![Hirai, Sato, Podsiadlowski, Vigna-Gomez, Mandel, 2020, MNRAS, 499, 1154](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.1154H/abstract)

Fixed mesh refinement <br>
![Hirai, Podsiadlowski, Owocki, Schneider, Smith, 2021, MNRAS, 503, 4276](https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.4276H/abstract)

Sink particles <br>
![Hirai & Podsiadlowski, 2022, MNRAS, 517, 4544](https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.4544H/abstract)
