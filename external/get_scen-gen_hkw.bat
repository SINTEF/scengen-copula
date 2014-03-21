@echo off

setlocal

set tgDir=scen-gen_HKW
set repo=https://kaut@bitbucket.org/kaut/code_scen-gen_hkw
set cbBin="c:\Program Files (x86)\CodeBlocks\codeblocks.exe"

if not exist %tgDir% (
	mkdir %tgDir%
	hg clone %repo% %tgDir%
)

if exist %cbBin% (
	%cbBin% --build --target=Release_dll --no-batch-window-close %tgDir%/scen-gen_HKW.cbp
	pause
	%cbBin% --build --target=Debug_dll --no-batch-window-close %tgDir%/scen-gen_HKW.cbp
)

endlocal