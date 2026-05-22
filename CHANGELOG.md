# Changelog

## [1.0.1] - 2026-05-26

* Regenerate slsqp.rs code from nlopt 2.10.1 slsqp C code with [c2rust 0.22.1](https://github.com/immunant/c2rust).
  Note: previous manual edits where re-introduced and libc dependency is removed

## [1.0.0] - 2026-01-14

* Bump version to 1.0.0
* Use Edition 2024

## [0.1.1] - 2023-10-19

* Fix default constraint tolerance (set to zero, was 2e-4) 

## [0.1.0] - 2023-10-13

SLSQP optimizer in Rust (SLSQP stands for Sequential Least SQuares Programming). 
SLSQP Rust code was generated from [NLopt](https://github.com/stevengj/nlopt) 2.7.1 SLSQP C code thanks 
to [c2rust](https://github.com/immunant/c2rust) transpiler then manually edited to make it work.

The algorithms is available as a `minimize` function. See documentation for further details.
Currently `slsqp` handle `xinit`, `bounds`, `max_eval`, `ftol_rel`, `ftol_abs`, `xtol_rel`, `xtol_abs` 
and only inequality contraints (like c <= 0). 
Gradients of the objective function and constraints have to be provided by the user for the method to work. 
