%module b1_func

%{
    #define SWIG_FILE_WITH_INIT
    #include "b1_func.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* x, int n)};

%include "b1_func.h"
