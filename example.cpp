#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"

using namespace alglib;
void  nlcfunc2_jac(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr)
{
    //
    // this callback calculates
    //
    //     f0(x0,x1,x2) = x0+x1
    //     f1(x0,x1,x2) = x2-exp(x0)
    //     f2(x0,x1,x2) = x0^2+x1^2-1
    //
    // and Jacobian matrix J = [dfi/dxj]
    //
    fi[0] = x[0]+x[1];
    fi[1] = x[2]-exp(x[0]);
    fi[2] = x[0]*x[0] + x[1]*x[1] - 1.0;
    jac[0][0] = 1.0;
    jac[0][1] = 1.0;
    jac[0][2] = 0.0;
    jac[1][0] = -exp(x[0]);
    jac[1][1] = 0.0;
    jac[1][2] = 1.0;
    jac[2][0] = 2*x[0];
    jac[2][1] = 2*x[1];
    jac[2][2] = 0.0;
}

int main(int argc, char **argv)
{
    //
    // This example demonstrates minimization of
    //
    //     f(x0,x1) = x0+x1
    //
    // subject to nonlinear inequality constraint
    //
    //    x0^2 + x1^2 - 1 <= 0
    //
    // and nonlinear equality constraint
    //
    //    x2-exp(x0) = 0
    //
    real_1d_array x0 = "[0,0,0]";
    real_1d_array s = "[1,1,1]";
    double epsx = 0.000001;
    ae_int_t maxits = 0;
    minnlcstate state;
    minnlcreport rep;
    real_1d_array x1;

    //
    // Create optimizer object and tune its settings:
    // * epsx=0.000001  stopping condition for inner iterations
    // * s=[1,1]        all variables have unit scale
    // * upper limit on step length is specified (to avoid probing locations where exp() is large)
    //
    minnlccreate(3, x0, state);
    minnlcsetcond(state, epsx, maxits);
    minnlcsetscale(state, s);
    minnlcsetstpmax(state, 10.0);

    //
    // Choose one of the nonlinear programming solvers supported by minnlc
    // optimizer:
    // * SLP - successive linear programming NLP solver
    // * AUL - augmented Lagrangian NLP solver
    //
    // Different solvers have different properties:
    // * SLP is the most robust solver provided by ALGLIB: it can solve both
    //   convex and nonconvex optimization problems, it respects box and
    //   linear constraints (after you find feasible point it won't move away
    //   from the feasible area) and tries to respect nonlinear constraints
    //   as much as possible. It also usually needs less function evaluations
    //   to converge than AUL.
    //   However, it solves LP subproblems at each iterations which adds
    //   significant overhead to its running time. Sometimes it can be as much
    //   as 7x times slower than AUL.
    // * AUL solver is less robust than SLP - it can violate box and linear
    //   constraints at any moment, and it is intended for convex optimization
    //   problems (although in many cases it can deal with nonconvex ones too).
    //   Also, unlike SLP it needs some tuning (penalty factor and number of
    //   outer iterations).
    //   However, it is often much faster than the current version of SLP.
    //
    // In the code below we set solver to be AUL but then override it with SLP,
    // so the effective choice is to use SLP. We recommend you to use SLP at
    // least for early prototyping stages.
    //
    // You can comment out line with SLP if you want to solve your problem with
    // AUL solver.
    //
    double rho = 1000.0;
    ae_int_t outerits = 5;
    minnlcsetalgoaul(state, rho, outerits);
    minnlcsetalgoslp(state);

    //
    // Set constraints:
    //
    // Nonlinear constraints are tricky - you can not "pack" general
    // nonlinear function into double precision array. That's why
    // minnlcsetnlc() does not accept constraints itself - only constraint
    // counts are passed: first parameter is number of equality constraints,
    // second one is number of inequality constraints.
    //
    // As for constraining functions - these functions are passed as part
    // of problem Jacobian (see below).
    //
    // NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
    //       linear and general nonlinear constraints. This example does not
    //       show how to work with boundary or general linear constraints, but you
    //       can easily find it in documentation on minnlcsetbc() and
    //       minnlcsetlc() functions.
    //
    minnlcsetnlc(state, 1, 1);

    //
    // Activate OptGuard integrity checking.
    //
    // OptGuard monitor helps to catch common coding and problem statement
    // issues, like:
    // * discontinuity of the target/constraints (C0 continuity violation)
    // * nonsmoothness of the target/constraints (C1 continuity violation)
    // * erroneous analytic Jacobian, i.e. one inconsistent with actual
    //   change in the target/constraints
    //
    // OptGuard is essential for early prototyping stages because such
    // problems often result in premature termination of the optimizer
    // which is really hard to distinguish from the correct termination.
    //
    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
    //            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
    //
    //            Other OptGuard checks add moderate overhead, but anyway
    //            it is better to turn them off when they are not needed.
    //
    minnlcoptguardsmoothness(state);
    minnlcoptguardgradient(state, 0.001);

    //
    // Optimize and test results.
    //
    // Optimizer object accepts vector function and its Jacobian, with first
    // component (Jacobian row) being target function, and next components
    // (Jacobian rows) being nonlinear equality and inequality constraints.
    //
    // So, our vector function has form
    //
    //     {f0,f1,f2} = { x0+x1 , x2-exp(x0) , x0^2+x1^2-1 }
    //
    // with Jacobian
    //
    //         [  +1      +1       0 ]
    //     J = [-exp(x0)  0        1 ]
    //         [ 2*x0    2*x1      0 ]
    //
    // with f0 being target function, f1 being equality constraint "f1=0",
    // f2 being inequality constraint "f2<=0". Number of equality/inequality
    // constraints is specified by minnlcsetnlc(), with equality ones always
    // being first, inequality ones being last.
    //
    alglib::minnlcoptimize(state, nlcfunc2_jac);
    minnlcresults(state, x1, rep);
    printf("%s\n", x1.tostring(2).c_str()); // EXPECTED: [-0.70710,-0.70710,0.49306]

    //
    // Check that OptGuard did not report errors
    //
    // NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
    //       1.0 to some of its components.
    //
    optguardreport ogrep;
    minnlcoptguardresults(state, ogrep);
    printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
    printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
    printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
    return 0;
}