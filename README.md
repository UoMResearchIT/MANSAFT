# MANSAFT
- [src](src) contains source code which uses open-source code from [MINPACK](http://www.netlib.org/minpack/)
   and [NL_OPT](https://nlopt.readthedocs.io/en/latest/) instead of the NAG library
- Execute `make` in the `src` directory to see an example of how to compile and run the software with gfortran
- The NL_OPT library must first be installed. See the [installation script](nl_opt/install-nlopt.sh) for details.
  `nl_opt/install-nlopt.sh` was based on 
  the official [installation instructions](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/)
  and can be run like this `bash install-nlopt.sh install_directory`
- To enable the compiler to find the NL_OPT library, the environment variable `LD_LIBRARY_PATH` needs to include
  the path to the `lib` directory of the NL_OPT installation. See [nl_opt/add-lib-path.sh](nl_opt/add-lib-path.sh) for an example.
