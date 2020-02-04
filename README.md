# MANSAFT
- Old code lives in the `old_code` directory. This calculates slowly and unreliably.
- The `NAG` directory uses functions from the NAG library to speed up calculations
- `minpack` is a rewrite of the NAG code using open-source code from [MINPACK](http://www.netlib.org/minpack/)
  to replace calls to the NAG library
  - Execute `make` in the `minpack` directory to see an example of how to compile and run the software with gfortran
