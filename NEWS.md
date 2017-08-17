VERSION - 0.4.1 - 07-04-17
- Unified documentation for all intrinsic models at a single source.
- TreEvo is now Windows-capable, by removing calls to quartz() and also allowing package doParallel to be used in replacement of package doMC, while potentially retaining multithreading functionality
	- doMC moved to Suggests, and doParallel added to Suggests
- NAMESPACE file converted over to now be automatically handled by roxygen2
- Added NEWS file for tracking changes to TreEvo
-Removed an apparent duplicate data file (simData), leaving simRun; modified references to simData accordingly
