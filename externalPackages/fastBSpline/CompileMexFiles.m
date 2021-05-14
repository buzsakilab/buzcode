function [] = CompileMexFiles()
    %Compiles mex files that accelerate B-spline evaluation
    mex evalBin.c
    mex evalBSpline.c
    mex evalBinTimesY.c
end