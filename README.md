# CarlemanLinearization.jl

| **Documentation** | **Status** | **Community** | **License** |
|:-----------------:|:----------:|:-------------:|:-----------:|
| [![docs-dev][dev-img]][dev-url] | [![CI][ci-img]][ci-url] [![codecov][cov-img]][cov-url] [![PkgEval][pkgeval-img]][pkgeval-url] [![aqua][aqua-img]][aqua-url] [![dev-commits][dev-commits-url]][dev-commits-target] | [![zulip][chat-img]][chat-url] [![JuliaHub][juliahub-img]][juliahub-url] | [![license][lic-img]][lic-url] |

[dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[dev-url]: https://juliareach.github.io/CarlemanLinearization.jl/dev/
[ci-img]: https://github.com/JuliaReach/CarlemanLinearization.jl/actions/workflows/test-master.yml/badge.svg
[ci-url]: https://github.com/JuliaReach/CarlemanLinearization.jl/actions/workflows/test-master.yml
[cov-img]: https://codecov.io/github/JuliaReach/CarlemanLinearization.jl/coverage.svg
[cov-url]: https://app.codecov.io/github/JuliaReach/CarlemanLinearization.jl
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/CarlemanLinearization.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/C/CarlemanLinearization.html
[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl
[dev-commits-url]: https://img.shields.io/github/commits-since/JuliaReach/CarlemanLinearization.jl/latest.svg
[dev-commits-target]: https://github.com/JuliaReach/CarlemanLinearization.jl
[chat-img]: https://img.shields.io/badge/zulip-join_chat-brightgreen.svg
[chat-url]: https://julialang.zulipchat.com/#narrow/stream/278609-juliareach
[juliahub-img]: https://juliahub.com/docs/General/CarlemanLinearization/stable/version.svg
[juliahub-url]: https://juliahub.com/ui/Packages/General/CarlemanLinearization
[lic-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[lic-url]: https://github.com/JuliaReach/CarlemanLinearization.jl/blob/master/LICENSE

This package implements the Carleman linearization transformation of polynomial
differential equations in Julia.

## Features

The following methods are available:

- Construction of the Carleman embedding using sparse matrices
- Explicit error bounds [FP17]
- Improved error bounds for dissipative systems [L20]

## Related libraries

- [carlin](https://github.com/mforets/carlin): Python library with similar functionality.
- [CarlemanBurgers](https://github.com/hermankolden/CarlemanBurgers): MATLAB implementation of the classical Carleman solution of the viscous Burgers equation, used in https://arxiv.org/abs/2011.03185.

## References

- [[FP17] Forets, Marcelo, and Amaury Pouly. Explicit error bounds for carleman linearization. arXiv preprint arXiv:1711.02552 (2017).](https://arxiv.org/abs/1711.02552)

- [[L20] Liu JP, Kolden HØ, Krovi HK, Loureiro NF, Trivisa K, Childs AM. Efficient quantum algorithm for dissipative nonlinear differential equations. arXiv preprint arXiv:2011.03185. 2020 Nov 6.](https://arxiv.org/abs/2011.03185)
