function options = nppcaOptions;

% NPPCAOPTIONS Set up an options vector for noisy PPCA.

% NPPCA

options.tol = 1e-4;
options.stepChecks = 0; % Whether to check likelihood after each update.
options.display = 0; % whether or not to display progress in a figure.
options.maxIters = 100; % Number of EM iterations.
options.optimiser = foptions; % optimisation options for Sigma.
options.optimiser(2) = options.tol;
options.optimiser(3) = options.tol;
