# slsqp - a pure Rust implementation

[![tests](https://github.com/relf/slsqp/actions/workflows/tests.yml/badge.svg)](https://github.com/relf/slsqp/actions/workflows/tests.yml)
[![crates.io](https://img.shields.io/crates/v/slsqp)](https://crates.io/crates/slsqp)
[![docs](https://docs.rs/slsqp/badge.svg)](https://docs.rs/slsqp)

SLSQP is a sequential quadratic programming algorithm for nonlinearly constrained gradient-based optimization based on the implementation by Dieter Kraft and described in:

  > Dieter Kraft, "A software package for sequential quadratic programming", Technical Report DFVLR-FB 88-28, Institut für Dynamik der Flugsysteme, Oberpfaffenhofen, July 1988.
  > Dieter Kraft, "Algorithm 733: TOMP–Fortran modules for optimal control calculations," ACM Transactions on Mathematical Software, vol. 20, no. 3, pp. 262-281 (1994).

The Rust code was generated/adapted from the C code from the [NLopt](https://github.com/stevengj/nlopt) project (version 2.7.1).
The algorithme is available here as a `minimize` function.
An initial transpilation was done with [c2rust](https://github.com/immunant/c2rust) then the code was manually edited to make it work. The callback mechanism is inspired from the Rust binding of NLopt, namely [rust-nlopt](https://github.com/adwhit/rust-nlopt).

## Example

```bash
cargo run --example paraboloid
```

## Related projects

* [rust-nlopt](https://github.com/adwhit/rust-nlopt): the Rust binding of the [NLopt project](https://nlopt.readthedocs.io)
* [cobyla](https://github.com/relf/cobyla): a pure Rust implementation of the COBYLA algorithm.

## License

The project is released under MIT License.
