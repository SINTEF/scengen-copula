@echo off

setlocal enabledelayedexpansion

set nFailed=0
set LOG_F=".\test_main.log"

:: ----------------------------------------------------------------------------
echo Testing the main generator
echo ==========================
set SG_BIN=".\scen-gen_cop"
set PARAMS_MAIN=(^
	"--cop-type normal --input normal_3d_corrs.dat --output normal_3d_cop.out 100"^
	"-c normal -i normal_3d_corrs.dat -o normal_3d_cop.out 100"^
	"--cop-type normal --input normal_3d_corrs.dat --marg-type normal --marg-par normal_3d_moms.dat --output normal_3d_val.out 100"^
	"-c normal -i normal_3d_corrs.dat -m normal -d normal_3d_moms.dat -o normal_3d_val.out 100"^
	"--cop-type normal --input normal_3d_corrs.dat --marg-type normal --marg-par normal_3d_moms.dat --output normal_3d_val.out --cop-output normal_3d_cop.out 100"^
	"-c normal -i normal_3d_corrs.dat -m normal -d normal_3d_moms.dat -o normal_3d_val.out --cop-output normal_3d_cop.out 100"^
	"--cop-type sample --input hist-data.dat --output test_cop.out 100"^
	"-c sample -i hist-data.dat -o test_cop.out 100"^
	"--cop-type sample --input hist-data.dat --marg-type sample --output test_val.out 100"^
	"-c sample -i hist-data.dat -m sample -o test_val.out 100"^
	"--cop-type sample --input hist-data.dat --marg-type normal --marg-par normal_3d_moms.dat --output test-normal_val.out 100"^
	"-c sample -i hist-data.dat -m normal -d normal_3d_moms.dat -o test-normal_val.out 100"^
	"--cop-type indep --dim 3 --marg-type moments --marg-par tg_moms.dat --mat-file-fmt 100"^
	"-c indep --dim 3 -m moments -d tg_moms.dat --mat-file-fmt 100"^
	"--cop-type normal --input normal_3d_corrs.dat --marg-type moments --marg-par tg_moms.dat --mat-file-fmt 100"^
	"-c normal -i normal_3d_corrs.dat -m moments -d tg_moms.dat --mat-file-fmt 100"^
	"--cop-type indep --dim 18 --marg-type mixed --marg-par mixed-margins.dat 100"^
	"-c indep --dim 18 -m mixed -d mixed-margins.dat 100"^
)
for %%p in %PARAMS_MAIN% do (
	set par=%%p
	set par=!par:"=!
	REM "
	echo TEST: !par!
	!SG_BIN! !par! > !LOG_F! 2>&1
	if ERRORLEVEL 1 (
		echo     : FAILED .. output follows:
		type !LOG_F!
		echo.
		set /a nFailed+=1
	) else (
		echo     : OK
	)
)

:: ----------------------------------------------------------------------------
echo.
echo Testing the forecast-error-based generator
echo ==========================================
set SG_BIN=".\scen-gen_cop_fc-err"
set PARAMS_MAIN=(^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --scens 6"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --periods 3 --scens 6"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --cur-val 2.5,6 --scens 6"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --fcast-incl-cur --scens 6"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --branching 4,3,2,1"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --branching 4,3,2"^
	"--hist-data hist-fcasts.dat --forecast forecast.dat --cur-val 2.5,6 --charts-per-var scen_for_var_ --branching 4,3,2,1"^
)
for %%p in %PARAMS_MAIN% do (
	set par=%%p
	set par=!par:"=!
	REM "
	echo TEST: !par!
	!SG_BIN! !par! > !LOG_F! 2>&1
	if ERRORLEVEL 1 (
		echo     : FAILED .. output follows:
		type !LOG_F!
		echo.
		set /a nFailed+=1
	) else (
		echo     : OK
	)
)

:: ----------------------------------------------------------------------------
:: report results
echo.
if !nFailed! EQU 0 (
	echo ALL TESTS PASSED
) else (
	echo WARNING: FAILED !nFailed! TESTS!
)

:: cleaning
if exist !LOG_F! (
	del !LOG_F!
)
