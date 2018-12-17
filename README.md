# methExprsHELIX
Analyze association between DNA Methylation and Expression in HELIX subcohort

Structure of the project:

* data: Raw data from the project
* results: Intermediate and result files
  + preprocessFiles: Objects with methylation and Gene expression to used in the analyses.
  + MethComBatExpResidualsCellAdj: Adjust for covariates and cell counts
  + MethComBatExpResidualsNoCellAdj: Adjust for covariates but not cell counts
* src: Scripts used to analyze the data