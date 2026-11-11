# Mathematical references

The following books provide standard foundations for the algorithms used in `alggeom`.

- David A. Cox, John Little, and Donal O'Shea (2015), *Ideals, Varieties, and Algorithms*, 4th ed., Springer. [Publisher page](https://link.springer.com/book/10.1007/978-3-319-16721-3). Gröbner bases, elimination, ideals, varieties, and the algebraic-geometric interpretation of computation.
- Gert-Martin Greuel and Gerhard Pfister (2008), *A Singular Introduction to Commutative Algebra*, 2nd ed., Springer. [Publisher page](https://link.springer.com/book/10.1007/978-3-540-73542-7). Primary decomposition, local algebra, standard bases, and computational constructions.
- David Eisenbud (1995), *Commutative Algebra: With a View Toward Algebraic Geometry*, Springer. [Publisher page](https://link.springer.com/book/10.1007/978-1-4612-5350-1). Dimension, multiplicity, homological algebra, normalization, and module-theoretic methods.
- Shuhong Gao (2003), *Factoring Multivariate Polynomials via Partial Differential Equations*, Mathematics of Computation 72(242), 801–822. [DOI](https://doi.org/10.1090/S0025-5718-02-01428-X). Exact logarithmic-derivative absolute factorization in characteristic zero.
- Wolfgang Ruppert (1986), *Reduzibilität ebener Kurven*, Journal für die reine und angewandte Mathematik 369, 167–191. [DOI](https://doi.org/10.1515/crll.1986.369.167). Geometric reducibility of plane curves.
- Volker Weispfenning (1992), *Comprehensive Gröbner Bases*, Journal of Symbolic Computation 14(1), 1–29. [DOI](https://doi.org/10.1016/0747-7171(92)90023-W). Parametric Gröbner bases under specialization.
- Antonio Montes (2002), *A New Algorithm for Discussing Gröbner Bases with Parameters*, Journal of Symbolic Computation 33(2), 183–208. [DOI](https://doi.org/10.1006/jsco.2001.0504). Constructible parameter branches and specialized bases.
- Deepak Kapur, Yao Sun, and Dingkang Wang (2013), *An Efficient Algorithm for Computing a Comprehensive Gröbner System of a Parametric Polynomial System*, Journal of Symbolic Computation 49, 27–44. [DOI](https://doi.org/10.1016/j.jsc.2011.12.015). Efficient branch refinement for comprehensive systems.

## Primary decomposition and parametric Gröbner bases

- Patrizia Gianni, Barry Trager, and Gail Zacharias (1988), *Gröbner Bases and Primary Decomposition of Polynomial Ideals*, *Journal of Symbolic Computation* **6**(2–3), 149–167. [DOI](https://doi.org/10.1016/S0747-7171(88)80040-3). Primary decomposition through elimination and factorization over finitely generated fields.
- Antonio Montes and Michael Wibmer (2010), *Gröbner Bases for Polynomial Systems with Parameters*, *Journal of Symbolic Computation* **45**(12), 1391–1425. [DOI](https://doi.org/10.1016/j.jsc.2010.06.017). Locally closed parameter strata with specialization-valid Gröbner data.
- Antonio Montes (2018), *The Gröbner Cover*, Algorithms and Computation in Mathematics **27**, Springer. [DOI](https://doi.org/10.1007/978-3-030-03904-2). Gröbner-cover algorithms and constructible-set representations, including Singular's `grobcov.lib`.

## Gröbner performance, modular algorithms, and change of ordering

- Elizabeth A. Arnold (2003), *Modular algorithms for computing Gröbner bases*, *Journal of Symbolic Computation* **35**(4), 403–419. [DOI](https://doi.org/10.1016/S0747-7171(02)00140-2). Modular computation over rational coefficient fields, unlucky-prime criteria, lifting, and exact verification.
- Jean-Charles Faugère, Patrizia Gianni, Daniel Lazard, and Teo Mora (1993), *Efficient computation of zero-dimensional Gröbner bases by change of ordering*, *Journal of Symbolic Computation* **16**(4), 329–344. [DOI](https://doi.org/10.1006/jsco.1993.1051). The FGLM algorithm, which converts Gröbner bases of zero-dimensional ideals between monomial orders.
- Rüdiger Gebauer and H. Michael Möller (1988), *On an installation of Buchberger's algorithm*, *Journal of Symbolic Computation* **6**(2–3), 275–286. [DOI](https://doi.org/10.1016/S0747-7171(88)80048-8). Redundant S-pair elimination and efficient Buchberger criteria.
- Werner M. Seiler (2002), *Taylor and Lyubeznik Resolutions via Gröbner Bases*, *Journal of Symbolic Computation* **34**(6), 597–608. [DOI](https://doi.org/10.1006/jsco.2002.0573). Schreyer syzygies, module normal forms, and Buchberger's chain criterion.
- Jean-Charles Faugère (1999), *A new efficient algorithm for computing Gröbner bases (F4)*, *Journal of Pure and Applied Algebra* **139**(1–3), 61–88. [DOI](https://doi.org/10.1016/S0022-4049(99)00005-5). Symbolic preprocessing and sparse linear-algebra reductions.
- Christian Eder and Jean-Charles Faugère (2017), *A survey on signature-based algorithms for computing Gröbner bases*, *Journal of Symbolic Computation* **80**(3), 719–784. [DOI](https://doi.org/10.1016/j.jsc.2016.07.031). Signature-safe reductions, rewriting criteria, and the differences between certified syzygy-based pair pruning and fully signature-based algorithms.
