Brief detail of various R files within the directory:

- getInfluenceFunctions.R: R file that takes an orm object and outputs influences functions for beta (and alpha if needed), that are then used for estimating GR weights and calculating the two-phase variance estimator. \

- pred_quantile.R: Function that returns the predicted quantile given a covariate profile, GR-weighted CPM, and estimated two-phase variance-covariance matrix. Based on internal code from Frank\'92s rms package.\

- pred_ExProb.R: Function that returns the predicted exceedance probability (Pr Y>y) given a covariate profile, GR-weighted CPM, and estimated two-phase variance-covariance matrix. Based on internal code from Frank\'92s rms package.\
