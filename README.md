## choibc

Compute High Order Impedance Boundary Condition (HOIBC) coefficients by a constrained optimisation problem

### Quick Install

```shell
    git clone https://github.com/pirpyn/choibc
```

### Quick Build

```shell
    make FC=g++ -j
```

### Description

This code is a C++ implementation of the [Fortran](https://github.com/pirpyn/fhoibc) which was written during my Ph.D thesis to provides a library to compute HOIBC coefficients that garantees the uniqueness of the solution of the exterior Maxwell equations.
To do so, we solve a constrained optimisation problem using SQP method.

Theses coefficient takes curvatures into account, as the local infinite cylinder and the sphere are implemented along with the classic infinite plane.

### Dependencies

This program needs 

  * [LAPACK](http://www.netlib.org/lapack/). It is recommended to use your distribution package manager to get `liblapack-dev` but a shell script `get_lapack.sh` is provided.

### Building

A Makefile is provided. To print all makefile recipes and options, type
```shell
    make help
```
### Testing

A sample program is provided in the `src/main` folder. This program
  * reads a json filename from the command line and initialise a `hoibc::data_t` structure accordingly.
  * pass this struct the `hoibc::main` function.

### Compiler support

This code has been tested with and should run with at least

  * `g++ 9.3.0

Newer version should compile and run fine.

### References

  * P. Payen, "Modélisation précise d'un matériaux hétérogène en hyperfréquence par une condition d'impédance d'ordre élevé", Ph.D Université Paris 13, Dec. 2019, http://theses.fr/s197772
  *  S. Bochkanov, "ALGLIB 3.16.0 for C++" (www.alglib.net)

### License

  * The original sourcecode and the modifications are released under a [BSD-3-Clause license](https://github.com/pirpyn/choibc/blob/master/LICENSE).
