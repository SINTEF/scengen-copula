# probability allocation model for the copula generation tool
#
# this version minimizes both the total and the max error
# plus the deviations from the initial (equiprobable) state

param nVars;
param nScens;

set VARS := {0 .. nVars-1};
set SCENS := {0 .. nScens-1};

param rank{VARS, SCENS} in SCENS;
param tgCdf{SCENS};

param defProb default 1.0 / nScens;
param minProbFrac default 0.0;

var prob{SCENS} >= minProbFrac * defProb <= 1.0;
var devPos{SCENS} >= 0;
var devNeg{SCENS} >= 0;
var maxDev;
var errPos{SCENS} >= 0;
var errNeg{SCENS} >= 0;
var maxErr;

s.t. cop_cdf{s in SCENS}:
	sum {ss in SCENS: forall{i in VARS} rank[i, ss] <= rank[i, s]} prob[ss]
	+ errPos[s] -errNeg[s] = tgCdf[s];

s.t. prob_def{s in SCENS}:
	prob[s] = defProb + devPos[s] - devNeg[s];

s.t. sum_prob:
	sum {s in SCENS} prob[s] == 1.0;

s.t. maxDev_err {s in SCENS}:
	devPos[s] + devNeg[s] <= maxDev;

s.t. maxErr_def {s in SCENS}:
	errPos[s] + errNeg[s] <= maxErr;

param wAvgErr default 1.0;
param wMaxErr default 1.0;
param wAvgDev default 1.0;
param wMaxDev default 1.0;
minimize tot_avg_error:
	wAvgErr * sum{s in SCENS} (errPos[s] + errNeg[s]) / nScens
	+ wMaxErr * maxErr
	+ wAvgDev * sum{s in SCENS} (devPos[s] + devNeg[s]) / nScens
	+ wMaxDev * maxDev
;

param scOutFile symbolic default "prob-alloc.out";

# -----------------------------------------------------

# broken in glpk 4.4
#display {s in SCENS} cop_cdf[s];

solve;

# command line output
printf "\n\nResults:\n";
printf "avg_error = %g\n", sum{s in SCENS} (errPos[s] + errNeg[s]) / nScens;
printf "max_error = %g\n", max{s in SCENS} (errPos[s] + errNeg[s]);
printf "avg_dev   = %g\n", sum{s in SCENS} (devPos[s] + devNeg[s]) / nScens;
printf "max_dev   = %g\n", max{s in SCENS} (devPos[s] + devNeg[s]);
printf "\nscenario-probabilities:\n";
printf {s in SCENS} "sc.%2d: %.4f\n", s, prob[s];
printf "\n";
printf "min prob = %g\n", min{s in SCENS} prob[s];
printf "\n";

# file output
printf "# optimal scenario probabilities from prob-allow.mod\n" > scOutFile;
printf "# %d scenarios\n", nScens >> scOutFile;
printf {s in SCENS} "%g\n", prob[s] >> scOutFile;

end;
