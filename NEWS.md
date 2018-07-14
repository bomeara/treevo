Version - 0.9.0 - 07-13-18
-Implemented new highest density region methods, replacing former credible interval functions

Version - 0.6.0 - 03-19-18
-Extensive modularization of doRun_prc code
-Full parellization of PRC algorithm for doRun_prc

Version - 0.5.0 - 01-14-18
-Seemed reasonable to bump version to 0.5 given great change in infrastructure

Version - 0.4.6 - 01-13-18
-Stable build with geiger::fitContinuous replaced for phylolm::phylolm
	-Stable tests, stable examples including dont-run examples
-Much faster for many package functionalities

Version - 0.4.5 - 01-11-18
-Reduced code redundancy in doRun_prc
-Many bug fixes
-Added many tests
-Prepared code to transition from using geiger::fitContinuous to phylolm::phylolm

VERSION - 0.4.2 - 10-03-17
- Internalized many functions with unclear need-to-use for typical users, and that were not necessary as input creators for other exported functions
- Internalized functions getSimulationSplits and getTaxonDFWithPossibleExtinction, and modified doSimulation functions to simply accept phylogenies as input rather than these secondary products
- Merged functions with shared design (particularly shared parameter sets) into the same help files to reduce redundancy
- Improved documentation and made working examples for all (although many are dont-test examples)

VERSION - 0.4.1 - 07-04-17
- Unified documentation for all intrinsic models at a single source.
- TreEvo is now Windows-capable, by removing calls to quartz() and also allowing package doParallel to be used in replacement of package doMC, while potentially retaining multithreading functionality
	- doMC moved to Suggests, and doParallel added to Suggests
- NAMESPACE file converted over to now be automatically handled by roxygen2
- Added NEWS file for tracking changes to TreEvo
-Removed an apparent duplicate data file (simData), leaving simRun; modified references to simData accordingly




