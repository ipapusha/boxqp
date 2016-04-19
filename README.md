# box_qp

## Requirements
This solver requires the command `ldlsparse(..)` 
which is available from the LDL package, a part of 
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html).

Here is an installation procedure:
```
localhost$ wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.5.2.tar.gz
localhost$ tar xvzf SuiteSparse-4.5.2.tar.gz
```
then, within Matlab:
```
> cd SuiteSparse/LDL/MATLAB/
> ldl_install
```
This will make the command `ldlsparse(..)` available for the solver, enabling
factorization caching. Use `pathtool` or `addpath` to add LDL to the Matlab
path.
