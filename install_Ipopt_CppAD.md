## Intalling Ipopt and CppAD

### Dependencies

At this point in the curriculum students will have set up their SDC Term 2 environment and dependencies, with the exception of Ipopt, Fortran, and CppAD.  If you are setting up a fresh environment please refer to setup instructions starting [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/382ebfd6-1d55-4487-84a5-b6a5a4ba1e47).

### Installation Process

1.  Clone this repository and navigate to the cloned directory
2.  [Download](https://www.coin-or.org/download/source/Ipopt/) the appropriate version of Ipopt (3.12.7 or higher) from the link below.  You may also use wget or a similiar command to download the source from the command line (see Linux instructions).
3.  Follow the instructions for your environment

* [Ipopt](https://projects.coin-or.org/Ipopt)
  * **Mac:**
    ```
      brew tap homebrew/science
      brew install ipopt --with-openblas
    ```
  * **Linux:**
    * ```sudo apt-get install gfortran```
    *  ```apt-get install unzip```
    * ```wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.7.zip && unzip Ipopt-3.12.7.zip && rm Ipopt-3.12.7.zip```
    * Call `install_ipopt.sh` with the source directory as the first argument, ex: ```./install_ipopt.sh Ipopt-3.12.7``` or ```bash install_ipopt.sh Ipopt-3.12.7```

  * **Windows:** For Windows environments there are two main options
    * Follow Linux instructions in the Ubuntu Bash environment
    * Use the docker container described [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77), which comes pre-configured with Ipopt.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * **Windows:** For Windows environments there are two main options
    * Follow Linux instructions in the Ubuntu Bash environment
    * Use the docker container described [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/16cf4a78-4fc7-49e1-8621-3450ca938b77), which comes pre-configured with CppAD.

### Mind the Line Quiz Dependencies
It may be neccesary to install additional dependencies, especially for Docker and Ubuntu BASH on Windows hosts.  A complete list of dependencies can be found [here](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile).  The mind the line solution uses a matplotlib inspired plotting cpp plotting library, which depends on python components.  To enable plotting, most Windows users will need to execute the following commands **note for dockers users, leave out ```sudo```):
- ```sudo apt-get update```
- ```sudo apt-get install python-matplotlib```
- ```sudo apt-get install python2.7-dev```

In addition, to display plots, an X-server must be running on the host and accessible.  To accomplish this in Ubuntu BASH for windows, do the following:
- [Download and install Xming](https://sourceforge.net/projects/xming/?source=typ_redirect)
- Start Xming
- execute the following in the terminal: ```export DISPLAY=:0```
- run the code from the build folder ```./mpc```

**Note to Docker Users** The Xming solution does not work out of the box since the docker container communicates with the VM and not the host.  To run the code without error it is necessary to comment out or remove the plotting code and the bottom of the MPC.cpp solution file.  A current work-around to visualizing results is to send the results to a file, transfer the file to the host, then use a visualaztion tool on the host.

### Troubleshooting

* If challenges to installation are encountered (install script fails).  Please consult the forums.  Please feel free to submit additional tips or forum threads to the [issue reports repo](https://github.com/udacity/sdc-issue-reports), for potential inclusion in this document.
*  **Some Mac users have experienced the following error:**
     ```
     Listening to port 4567
     Connected!!!
     mpc(4561,0x7ffff1eed3c0) malloc: *** error for object 0x7f911e007600: incorrect checksum for freed object
     - object was probably modified after being freed.
     *** set a breakpoint in malloc_error_break to debug
     ```
     This error has been resolved by updrading ipopt with
     ```brew upgrade ipopt --with-openblas```
     per this [forum post](https://discussions.udacity.com/t/incorrect-checksum-for-freed-object/313433/19)

** Useful Resources:**

- [A reference url for general GUI use from Ubuntu BASH](https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/)
- [A helpful troubleshooting thread](https://discussions.udacity.com/t/what-call-to-subplot-failed/298481)
- [A Mind the Line issue and solution collection thread](https://discussions.udacity.com/t/error-loading-module-matplotlib-pyplot/248132/8?source_topic_id=298481)
- [Matplotlibcpp](https://github.com/lava/matplotlib-cpp/tree/master/examples)
