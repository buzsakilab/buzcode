function runningOctave = isOctave()

% Returns whether this code has been run using Octave or MATLAB


try 
    OCTAVE_VERSION;
    runningOctave = 1;
catch
    runningOctave = 0;    
end














