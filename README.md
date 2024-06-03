# Rank-Revealing LU and QR Factorizations

This is the code repository accompanied the manuscript titled "[How to reveal the rank of a matrix?](https://arxiv.org/abs/2405.04330)" It implements four major functions:

* [A rank-revealing (partial) LU factorization](./RRLU.m)
* [A rank-revealing (partial) QR factorization](./RRQR.m)
* [A volume ratio metric computation for LU factorizations](./gammaLU.m)
* [A volume ratio metric computation for QR factorizations](./gammaQR.m)

All four functions are well-documented and can be queried using the MATLAB `help` function. In addition, this repository contains two scripts [Example_LU.m](./Example_LU.m) and [Example_QR.m](./Example_QR.m) that contain examples of function calls. The code was implemented and tested using MATLAB 2024b.