# SSAF_2024_repo

<!--- Created on 2024.06.24.
Last change on 2024.06.24. -->

This is a collection of MATLAB files associated with the preprint [A single shooting method with approximate Fréchet derivative for computing geodesics on the Stiefel manifold](https://arxiv.org/abs/2404.04089).

In that paper, we propose a single shooting method with approximate Fréchet derivative (SSAF) to compute the Riemannian distance on the Stiefel manifold St(n,p), the set of n-by-p
matrices with orthonormal columns. The proposed method is a shooting method in the sense of the classical shooting methods for solving boundary value problems; see, e.g., Stoer and Bulirsch, 1991. The main feature is that we provide an approximate formula for the Fréchet derivative of the geodesic involved in our shooting method.

If there are any problems or bugs, feel free to [contact me](mailto:msutti@ncts.ntu.edu.tw).

If you use SSAF for your research, please kindly cite the following preprint:

> M. Sutti, 
"[A single shooting method with approximate Fréchet derivative for computing geodesics on the Stiefel manifold](https://arxiv.org/abs/2404.04089)." Tech. report, April 2024, https://arxiv.org/abs/2404.04089.

## I) Version History

- Ver 1, 25 Jun 2024: initial release.

## II) Contents

The main folder `SSAF_2024_repo` contains the following subfolders:

- `plots`: the folder where the plots are saved.
- `postprocessing`: contains the utilities for plotting and visualizing data.
- `utilities`: contains all the utility functions shared by all drivers.

## III) Installation and Usage

No installation is required. Use the `Driver_*.m` files to run the calculations or generate the figures appearing in the preprint.
The figures will then be saved into the `plots` folder.

Here is a list of the `Driver_*.m` files:

- `Driver_Bryner_method.m` runs the method of Bryner from the paper [Endpoint Geodesics on the Stiefel Manifold Embedded in Euclidean Space](https://epubs.siam.org/doi/10.1137/16M1103099), SIAM J. Matrix Anal. Appl., Vol. 38 (2017), No. 4, pp. 1139–1159;
- `Driver_comparison_methods.m` performs comparisons with existing state-of-the-art methods;
- `Driver_expm_skewsymm_matrix.m` calculates the averaged time for computing a matrix exponential expm of a unit norm 1000×1000 skew-symmetric matrix;
- `Driver_Fig4_Bryner.m` tries to reproduce Fig. 4 from Bryner's paper cited above;
- `Driver_Shooting_Approx_Frechet.m` runs an instance of our SSAF method.


## IV) License

The code written in this repository is GPL licensed.

## V) Dependencies

The file `Stiefel_Log_supp.m` in the `utilities` folder is taken from Appendix C of the paper [A Matrix-Algebraic Algorithm for the Riemannian Logarithm on the Stiefel Manifold under the Canonical Metric](https://epubs.siam.org/doi/10.1137/16M1074485), SIAM J. Matrix Anal. Appl., Vol. 38 (2017), No. 2, pp. 322–342.

