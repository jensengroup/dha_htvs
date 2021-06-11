# High-Throughput Virtual Screening of DHA/VHF
Code for the article: [High throughput virtual screening of 230 billion molecular solar heat battery candidates](https://doi.org/10.7717/peerj-pchem.16). The code is divided into four directories, each performing one task.

1. [run_sqm_calculations](run_sqm_calculations): contains code to compute SQM barrier heights and storage densities automatically.
2. [train_ml_models](train_ml_models): includes two notebooks to train ML models:
    * Train linear regression to predict storage densities for all 230 billion molecules.
    * Train LightGBM models to predict storage densities/barriers for 420 million molecules.
3. [run_all_230b](run_all_230b): is the Rust code used to run predictions on all 230 billion molecules, using the [linear regression model](train_ml_models/train_linear_ML_models.ipynb).
4. [run_ga_test](run_ga_test): is a notebook used to run the Genetic Algorithm (GA) test.
