Copula-based scenario generation
================================

The code implements copula-based heuristic from paper
'<em>A copula-based heuristic for scenario generation</em>' by Michal Kaut,
published in <em>Computational Management Science</em>, 11 (4) pp. 503-516, 2014,
[doi:10.1007/s10287-013-0184-4](http://dx.doi.org/doi:10.1007/s10287-013-0184-4). 


Basic usage
-----------

### Specifying the copula

The copula type/family is specified by option \c \--cop-type; its parameters are then read from file given by option \c \--input, with default value `cop-params.dat'.

At the moment, the code handles the following copula types:

- normal - the param. file contains the correlation matrix
- sample - the param. file contains a matrix of sample/historical values
            that the scenarios should replicate
- indep - independent random variables; the param. file contains only the
            dimension - an alternative is to give this directly using option
            \c \--dim, then the input file won't be read

### Specifying marginal distributions

If we do not specify any marginal distribution, the code will output
the copula itself, possibly scaled to ranks using option \c \--cop-as-ranks.

TO DO


License
-------

The code is dual-licensed under the
[GNU Lesser General Public License (LGPL), version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html) and the [Eclipse Public License (EPL)](http://www.eclipse.org/legal/epl-v10.htm).

For information about the licenses, including compatibility with other licenses,
see the official [EPL FAQ](http://www.eclipse.org/legal/eplfaq.php),
the official [GNU licenses compatibility matrix](https://www.gnu.org/licenses/gpl-faq.html#AllCompatibility),
or the relevant Wikipedia pages: 
[EPL](http://en.wikipedia.org/wiki/Eclipse_Public_License") and 
[LGPL](https://en.wikipedia.org/wiki/GNU_Lesser_General_Public_License).

Note that the EPL is not GPL-compatible. On the other hand, it is the license used in most [COIN-OR](https://www.coin-or.org/) projects.

One important aspect of both licenses (so called "weak copyleft") is that if you make any modification or addition to the code itself, you _must_ put your modification under the same license, the LGPL or EPL.

Note that it is explicitly <b>not needed</b> to put any application under the LGPL or EPL, if that application is just _using_ the code, without making any changes to it.

Parts of the code may be available under other licenses that are compatible with the LGPL and the EPL. For instance, source code that is governed by licenses like BSD or MIT, or source code that has been put into the Public Domain by its author can be used within the code. Files or folders including such source code mention those other licenses explicitly.
