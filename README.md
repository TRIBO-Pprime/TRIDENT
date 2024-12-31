Brief description
-----------------

**TRIDENT** - **TRI**bologie **DENT**aire - is designed to characterize dental surface microwear. It:

* reads a batch file
* transforms each keyword in an action
* when encountered, performs an analysis
* and writes the results in a specific file

The programs are written in Fortran (2003+)

[top](#table-of-contents)

Building TRIDENT
--------------

Make sure that the dependencies described below are built.

Debug mode:

```bash
make debug
```
Normal mode:

```bash
make
```

Dependencies
------------

TRIDENT needs some components that are available through the following packages:

* *TOOLIB* ... Some general tools like FFT, file handling, minimization, sorting, etc.

* *TPGLIB* ... Some more specific programs like filtering, anisotropy analysis, derivation, etc.

As a the packages [TOOLIB](https://github.com/TRIBO-Pprime/TOOLIB) and [TPGLIB](https://github.com/TRIBO-Pprime/TPGLIB) have to be downloaded.

[top](#table-of-contents)

Third party components
----------------------

TOOLIB also uses external codes such as:

+ *FFTW3.3* ... [Fastest Fourier Transform in the West](https://www.fftw.org/) ... GNU General Public License

+ *Pikaia_oop* ... [Modern Fortran Edition of the Pikaia Genetic Algorithm](http://github.com/jacobwilliams/pikaia) ... BSD like

+ *GNUFOR* ... [Gnuplot Fortran](https://people.math.sc.edu/Burkardt/f_src/gnufor/gnufor.html) ... GNU General Public License

+ *Bspline-fortran* ... [Multidimensional B-Spline Interpolation of Data on a Regular Grid](https://github.com/jacobwilliams/bspline-fortran) ... BSD like

`.sur` surface files can be visualized and analyzed with [Gwyddion software](http://gwyddion.net/download.php), a modular program for SPM (scanning probe microscopy) data visualization and analysis.

Typical use
-----------

The program `main` reads a script file `my_script.f90` where the following parameters are typically defined:

+ the number of threads
+ the surface to be analyzed (.sur file)
+ the sampling parameters
+ the analyses to be performed

Run:

```bash
./main cfg/my_script.md
```

[top](#table-of-contents)

Example of use
--------------

## Test 1 -- 1 surface, global analysis

Run:

```bash
./main cfg/test1.job
```

## Test 2 -- 2 surfaces, global and sampling analysis

Run:

```bash
./main cfg/test2.job
```

Documentation
-------------
The documentation is automatically generated with [FORD](https://github.com/Fortran-FOSS-Programmers/ford), an automatic documentation generator for modern Fortran programs.

[top](#table-of-contents)

License
-------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warrenty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

You should have received a [copy](https://github.com/TRIBO-Pprime/TRIDENT/LICENSE) of the GNU General Public License along with this program. If not, see the [GNU website](https://www.gnu.org/licenses/gpl.html).

