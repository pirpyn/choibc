## choibc

Compute High Order Impedance Boundary Condition (HOIBC) coefficients by a constrained optimisation problem

### Quick Install

```shell
git clone https://github.com/pirpyn/choibc
```

### Quick Build

```shell
make CXXC=g++ -j
```

### Description

This code is a C++ implementation of the [Fortran](https://github.com/pirpyn/fhoibc) which was written during my Ph.D thesis to provides a library to compute HOIBC coefficients that garantees the uniqueness of the solution of the exterior Maxwell equations.
To do so, we solve a constrained optimisation problem using SQP method.

Theses coefficient takes curvatures into account, as the local infinite cylinder and the sphere are implemented along with the classic infinite plane.

### Dependencies

This program needs

  * [LAPACKE](http://www.netlib.org/lapack/). On Ubuntu, the apt command is `apt get liblapacke-dev`.

### Building

A Makefile is provided. To print all makefile recipes and options, run
```shell
make help
```
### Testing

A `src/test` folder is present and check accuracy of bessel functions, the overloaded `operator*` and the `hoibc::norm` function.
To launch thoses test, run
```shell
make run_test
```

A sample program is provided in the `src/main` folder. This program
  * reads a json filename from the command line and initialise a `hoibc::data_t` structure accordingly,
  * pass this struct the `hoibc::main` function,
  * compute the errors between the IBC and the exact impedance operator.

To run it, launch
```shell
make run ARGS=input.json
```

### Compiler support

This code has been tested with and should run with at least

  * `g++ 9.3.0`
  * `clang 10.0`

Newer version should compile and run fine.

### References
  * P. Payen, O. Lafitte, B Stupfel, "A high order impedance boundary condition with uniqueness conditions to solve thetime harmonic Maxwell’s equations", Proceedings 	14th International Conference on Mathematical and Numerical Aspects of Wave Propagation. Book of Abstracts, https://doi.org/10.34726/waves2019, [direct link](  https://repositum.tuwien.at/bitstream/20.500.12708/637/2/th%20International%20Conference%20on%20Mathematical%20and%20Numerical%20Aspects%20of%20Wave%20Propagation%20Book%20of%20Abstracts.pdf#157_abstract)
  * P. Payen, "Modélisation précise d'un matériaux hétérogène en hyperfréquence par une condition d'impédance d'ordre élevé", Ph.D Université Paris 13, Dec. 2019, http://theses.fr/s197772
  *  S. Bochkanov, "ALGLIB 3.16.0 for C++" (www.alglib.net)

### License

  * The original sourcecode and the modifications are released under a [MIT license](https://raw.githubusercontent.com/pirpyn/choibc/master/LICENCE).
