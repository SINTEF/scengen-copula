#!/bin/sh

tgDir=scen-gen_HKW
repo=https://kaut@bitbucket.org/kaut/code_scen-gen_hkw
cbBin=codeblocks

if [ ! -d ${tgDir} ]; then
	mkdir $tgDir
	hg clone ${repo} ${tgDir}
fi

# this prints the path to $cbBin and returns 0 if OK, >0 if not found
command -v $cbBin >/dev/null 2>&1
if [ $? -eq 0 ]; then
	$cbBin --build --target=Debug_dll --no-batch-window-close $tgDir/scen-gen_HKW.cbp
	$cbBin --build --target=Release_dll --no-batch-window-close $tgDir/scen-gen_HKW.cbp
fi
