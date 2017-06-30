The directory contains additional code to make the NDT work with Octave.

The directory contains a version of the union.m function that will allow one to take the union
of a cell array and an empty vector (this works in MATLAB but not in Octave).
The private directory contains a helper function to get the union.m function 
to work. 


Useful note: if you want to save results run in Octave and load them in MATLAB then use
the '-mat7-binary' flag, e.g., save('my_results.mat', '-mat7-binary', 'DECODING_RESULTS')









