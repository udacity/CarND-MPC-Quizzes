# CarND Controls Quizzes

Quizzes for *Vehicle Models* and *Model Predictive Control* sections.

1. [Global Kinematic Model](./global_kinematic_model) - Implement the *Global Kinematic Model*.
2. [Polynomial Fitting](./polyfit) - Fit and evaluate polynomials.
3. [MPC](./mpc_to_line) - Implement MPC and minimize cross track and orientation errors to a straight line trajectory.

To do a quiz:

1. Go to quiz directory.
2. Make a build directory with `mdkir build`.
3. Change into the build directory, `cd build`.
4. Compile the project, `cmake .. && make`.

A solution for each quiz is presented in the solution directory.

## Dependencies

The *Global Kinematic Quiz* and *Polynomial Fitting* quizzes have all the dependencies in repo. For the *MPC* quiz
you'll have to install Ipopt and CppAD.

* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt --with-openblas`
  * Linux
    * [This Dockerfile](./Dockerfile) might be helpful.
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x, sudo apt-get install coinor-libipopt-dev. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from [here](https://www.coin-or.org/download/source/Ipopt/) or [here](https://github.com/coin-or/Ipopt/releases/) if the former doesn't work.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
