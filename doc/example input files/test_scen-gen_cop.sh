#!/bin/bash

# array for all the tests
declare -a ALL_TESTS=()

# ---------------------------------------------------------------------
# Testing the main executable
SG_BIN="./scen-gen_cop"
declare -a PARAMS_MAIN=(
	# 1. Generate 100 scenarios from normal copula with correlations from file
	"--cop-type normal --input normal_3d_corrs.dat --output normal_3d_cop.out 100"
	"-c normal -i normal_3d_corrs.dat -o normal_3d_cop.out 100"
	# 2. Generate 100 scenarios from normal distribution with correlations from file normal_3-vars_corrs.dat and means and variances from normal_3-vars_moms.dat,
	#    with output to file normal_3d_val.out
	"--cop-type normal --input normal_3d_corrs.dat --marg-type normal --marg-par normal_3d_moms.dat --output normal_3d_val.out 100"
	"-c normal -i normal_3d_corrs.dat -m normal -d normal_3d_moms.dat -o normal_3d_val.out 100"
	# 3. Like 2., but output also the copula, to file normal_3d_cop.out
	"--cop-type normal --input normal_3d_corrs.dat --marg-type normal --marg-par normal_3d_moms.dat --output normal_3d_val.out --cop-output normal_3d_cop.out 100"
	"-c normal -i normal_3d_corrs.dat -m normal -d normal_3d_moms.dat -o normal_3d_val.out --cop-output normal_3d_cop.out 100"
	# 4. Generate 100 scenarios of copula from file hist-data.dat, to file test_cop.out
	"--cop-type sample --input hist-data.dat --output test_cop.out 100"
	"-c sample -i hist-data.dat -o test_cop.out 100"
	# 5. Like 4., but use also margins from hist-data.dat and output to test_val.out
	"--cop-type sample --input hist-data.dat --marg-type sample --output test_val.out 100"
	"-c sample -i hist-data.dat -m sample -o test_val.out 100"
	# 6. Generate 100 scenarios with copula from hist-data.dat and normal margins with moments from normal_3d_moms.dat, with output to test-normal_val.out
	"--cop-type sample --input hist-data.dat --marg-type normal --marg-par normal_3d_moms.dat --output test-normal_val.out 100"
	"-c sample -i hist-data.dat -m normal -d normal_3d_moms.dat -o test-normal_val.out 100"
	# 7. Generate 100 scenarios of two independent margins with moments from tg_moms.dat, with moments compatible with spreadsheet formulas
	"--cop-type indep --dim 3 --marg-type moments --marg-par tg_moms.dat --mat-file-fmt 100"
	"-c indep --dim 3 -m moments -d tg_moms.dat --mat-file-fmt 100"
	# 8. Simulating the old moment-matching heuristic (except that the correlations are approximate):
	#   generate 100 scenarios from normal copula with correlations from normal_3d_corrs.dat and moments from tg_moms.dat, using the matrix input format
	"--cop-type normal --input normal_3d_corrs.dat --marg-type moments --marg-par tg_moms.dat --mat-file-fmt 100"
	"-c normal -i normal_3d_corrs.dat -m moments -d tg_moms.dat --mat-file-fmt 100"
	# 9. Generating 100 scenarios from independent mixed margins, showcasing all available margin type
	"--cop-type indep --dim 18 --marg-type mixed --marg-par mixed-margins.dat 100"
	"-c indep --dim 18 -m mixed -d mixed-margins.dat 100"
	# X. Generating scenarios with normal copula and margins from a sample is not yet implemented!!!
	#"--cop-type normal --input normal_3d_corrs.dat -m sample -d hist-data.dat --output normal-test_val.out 100"
)
for par in "${PARAMS_MAIN[@]}"; do
	ALL_TESTS=("${ALL_TESTS[@]}" "$SG_BIN $par")
done


# ---------------------------------------------------------------------
# Testing the library driver (does not take parameters)
SG_BIN="./scen-gen_cop_ex"
ALL_TESTS=("${ALL_TESTS[@]}" "$SG_BIN")
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH


# ---------------------------------------------------------------------
# Testing the forecast-error-based code
SG_BIN="./scen-gen_cop_fc-err"
PARAMS_FCERR=(
	# 1. basic syntax for two-stage tree; no root
	"--hist-data hist-fcasts.dat --forecast forecast.dat --scens 6"
	# 2. like 1., but limiting the horizon to 3 periods
	"--hist-data hist-fcasts.dat --forecast forecast.dat --periods 3 --scens 6"
	# 3. two-stage tree; explicit current value
	"--hist-data hist-fcasts.dat --forecast forecast.dat --cur-val 2.5,6 --scens 6"
	# 4. two-stage tree; current value from the forecast matrix
	"--hist-data hist-fcasts.dat --forecast forecast.dat --fcast-incl-cur --scens 6"
	# multi-stage tree
	"--hist-data hist-fcasts.dat --forecast forecast.dat --branching 4,3,2,1"
	# multi-stage tree with limited horizon
	"--hist-data hist-fcasts.dat --forecast forecast.dat --branching 4,3,2"
	# multi-stage tree with specified current value and output figures (file scen_for_var_X.png)
	"--hist-data hist-fcasts.dat --forecast forecast.dat --cur-val 2.5,6 --charts-per-var scen_for_var_ --branching 4,3,2,1"
)
for par in "${PARAMS_FCERR[@]}"; do
	ALL_TESTS=("${ALL_TESTS[@]}" "$SG_BIN $par")
done

LOG_F="./test_main.log"
nFailed=0
for test in "${ALL_TESTS[@]}"; do
	echo "TEST: $test"
	$test &> $LOG_F
	if [ $? == 0 ]; then
		echo "    : OK"
	else
		echo "    : FAILED .. output follows:"
		cat $LOG_F
		echo
		let nFailed=nFailed+1
	fi
done

if [ $nFailed == 0 ]; then
	echo "ALL TESTS PASSED"
else
	echo "WARNING: FAILED $nFailed TESTS!"
fi

# cleaning
if [ -f $LOG_F ]; then
	rm $LOG_F
fi
