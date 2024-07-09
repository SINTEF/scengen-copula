Copula-based scenario generation
================================

The code implements copula-based heuristic from paper
'<em>A copula-based heuristic for scenario generation</em>' by Michal Kaut,
published in <em>Computational Management Science</em>, 11 (4) pp. 503-516, 2014,
[doi:10.1007/s10287-013-0184-4](http://dx.doi.org/doi:10.1007/s10287-013-0184-4). 


Basic usage - `scen-gen_cop` binary
-----------------------------------

This is the main binary.
For exact syntax, run it with `--help` and/or check the examples in the `test_scen-gen_cop` script in `doc/example input files` and/or the examples in `doc/examples of use.txt`.

### Specifying the copula

The copula type/family is specified by option `--cop-type`; its parameters are then read from file given by option `--input`, with default value `cop-params.dat`.

At the moment, the code handles the following copula types:

- `normal` - the param. file contains the correlation matrix
- `sample` - the param. file contains a matrix of sample/historical values
             that the scenarios should replicate
- `indep` - independent random variables; the param. file contains only the
            dimension - an alternative is to give this directly using option
            `--dim`, then the input file won't be read
- `mixed` - the param. file contains a list of copula specifications per variable pair,
            given as `i j COPULA_TYPE [COP_PARAMS]`.
    - for exact syntax, see `cop_mixed_4d.dat` in `doc/example input files`.

### Specifying marginal distributions

Margins are specified using the `--marg-type` option, with values taken from file given by `--marg-par`.
Currently, the margin type is one of:

- `normal` for normal distrition
    - the param. file includes the means and standard deviations of the margin
- `moments` for margins given by their first 4 moments
    - the param. file includes the first 4 moments of the margin, one on each line
- `sample` for marginal distributions taken from historical data
    - at the moments, this assumes that copula is also of type `sample` and the data will be read from the file specified by `--input`
- `mixed` for marginal distributions specified individually for each margin

If we do not specify any marginal distribution, the code will output
the copula itself, possibly scaled to ranks using option `--cop-as-ranks`.


Forecast-based generator - `scen-gen_cop_fc-err` binary
-------------------------------------------------------

This is a forecast-based generator described in [this unpublished report](https://work.michalkaut.net/#Kaut2016).

It builds a multi-stage scenario tree based on:

- historical data
- historical forecast errors
    - these may be computed from historical forecasts and data
- current forecast
- current value (for 'now'), if not included in the forecast
- branching structure of the tree to be generated

To see the exact syntax, run the code with `--help` and/or check the examples in the `test_scen-gen_cop` script in `doc/example input files`.


Building
--------

### Dependencies

The code uses several [Boost](https://www.boost.org/) C++ libraries, which must be available in the include path.
While most Boost libraries are header-only, the code uses the [Program Options](https://www.boost.org/doc/libs/1_85_0/doc/html/program_options.html) library that must be installed.
This can be done with: 

- in MSYS2 with UCRT64 env:
    - run `pacman -S mingw-w64-ucrt-x86_64-boost` from the UCRT64 terminal
    - this installs the complete boost, so it may be too much
- in Ubuntu
    - install system package `libboost-all-dev` (for the complete Boost libraries) or `libboost-program-options-dev` for only the program-options library. The latter should be enough, but was not tested.


### Using Code::Blocks

The repository includes a project file `copula-gen.cbp` for the [Code::Blocks](https://www.codeblocks.org/) IDE. This is an open-source IDE, available for both Windows and Linux.

The most relevant build targets are:

- `Release` that makes the main `scen-gen_cop` binary
- `fc-err_Rel` that makes the forecast-based `scen-gen_cop_fc-err` binary
- `lib_Release` that makes the `scen-gen_cop` shared library (dll/so).

The targets can also be built from the terminal using

```console
codeblocks copula-gen.cbp --build --target=<TARGET>
```

### Using makefiles

The repository includes three platform-specific makefiles for the GNU make tools, generated from the `copula-gen.cbp` project file using `cbp2make` (included in Code::Blocks).

The target names are sanitized in the makefiles, so the above targets become `release`, `fc_err_rel`, and `lib_release`, respectively.

For example, on Linux, the command to build the release target is

```console
make.exe -f Makefile.unit release
```

On Windows with MSYS2, this becomes

```console
mingw32-make.exe -f Makefile.win release
```

**NB:** the Apple-specific `Makefile.mac` has not been tested at all.

**NB:** all the makefiles use the same paths, so one needs to clean the build environment (esp. the `.obj` directory) between builds for multiple platforms!

### Using Visual Studio

The repository includes project files for Visual Studio, but they do not work at the moment.

(The NuGet package manager finds only version 1.63 of `boost_smart_ptr-src`, while the rest of the Boost libraries is at version 1.85...)


### Tested platforms

The code is currently being built with GCC on the following platforms

- [MSYS2](https://www.msys2.org/) tools with the (default) UCRT64 environment.
    - install from [https://www.msys2.org] or with `winget install MSYS2.MSYS2`
	- if the install fails, try to disable some Windows security as described in [https://github.com/appveyor/ci/issues/3777] and [https://github.com/msys2/msys2.github.io/issues/303] - the security can be re-activated after installation
    - install the following packages (using `pacman -S <PACKAGE>`):
        - `mingw-w64-ucrt-x86_64-toolchain`
        - `mingw-w64-ucrt-x86_64-boost`
    - **NB:** the resulting binaries must be able to access (the correct version of) `libstdc++-6.dll`, otherwise they 'die' without any warning the first time they use some function from the STD library!
        - this could be avoided using the `-static` option for the linker, at the cost of larger binaries
- Ubuntu on WSL
    - requires only installation of the boost libraries, as described above


License
-------

The code is licensed under the [Mozilla Public License, version 2.0](https://www.mozilla.org/en-US/MPL/2.0/).
For more information, see the [official FAQ](https://www.mozilla.org/en-US/MPL/2.0/FAQ/) or the [Wikipedia article](https://en.wikipedia.org/wiki/Mozilla_Public_License) about the license.
