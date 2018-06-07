# `isat_ffd`
`isat_ffd`  is a reduced order model trained by FFD simulations that can be used for fast predictions of airflow in the building. For the detailed performance evaluation of `isat_ffd`, one should refer to the paper titled [Fast and Self-Learning Indoor Airflow Simulation Based on In Situ Adaptive Tabulation](https://www.colorado.edu/lab/sbs/sites/default/files/attached-files/2018_wtian_isat_1.pdf).

`ffd`, which is the abbreviation of Fast Fluid Dynamics, is parallelized in OpenCL to run full-scale simulations on multi-core devices, and generates training data for the reduced order model `isat`. For the performance of `ffd` parallelized in OpenCL, one can refer to the paper titled [A Systematic Evaluation of Accelerating Indoor Airflow Simulations Using Cross Platform Parallel Computing](https://www.colorado.edu/lab/sbs/sites/default/files/attached-files/2016_wtian_ffd_gpu.pdf)

`isat` was originally developed by Professor Stephen B. Pope at the Cornell University to efficiently simulate the racting flow with detailed chemistry. For more information regarding `isat`, one can go to the [webpage](https://tcg.mae.cornell.edu/isat.html)

# how to compile
TBA
# how to run exes
TBA
# who do i talk to
Wei Tian, Wei.Tian@Schneider-Electric.com
