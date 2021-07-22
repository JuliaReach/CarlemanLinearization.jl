# CarlemanLinearization.jl

[![Build Status](https://github.com/JuliaReach/CarlemanLinearization.jl/workflows/CI/badge.svg)](https://github.com/JuliaReach/CarlemanLinearization.jl/actions?query=workflow%3ACI)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliareach.github.io/CarlemanLinearization.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/juliareach/CarlemanLinearization.jl/blob/master/LICENSE)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This package implements the Carleman linearization transformation of polynomial
differential equations in Julia.

## Features

The following methods are available:

- Construction of the Carleman embedding using sparse matrices
- Explicit error bounds [FP17]
- Improved error bounds for dissipative systems [L20]

## Example

```julia
using CarlemanLinearization

....
```

## Related libraries

- [carlin](https://github.com/mforets/carlin): Python library with similar functionality.
- [CarlemanBurgers](https://github.com/hermankolden/CarlemanBurgers): MATLAB implementation of the classical Carleman solution of the viscous Burgers equation, used in https://arxiv.org/abs/2011.03185.

## References

- [FP17] Forets, Marcelo, and Amaury Pouly. "[Explicit error bounds for carleman linearization.](https://arxiv.org/abs/1711.02552)" arXiv preprint arXiv:1711.02552 (2017).

- [L20] Liu JP, Kolden HØ, Krovi HK, Loureiro NF, Trivisa K, Childs AM. "[Efficient quantum algorithm for dissipative nonlinear differential equations.](https://arxiv.org/abs/2011.03185). arXiv preprint arXiv:2011.03185. 2020 Nov 6.
