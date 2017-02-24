function logger(string, nlog)
  
  global LOGGER_LOGS

  
  print_on_stdout = 1;
  
  if nargin < 2 
    nlog = 1;
  end
  
  str = [datestr(clock, 31) ' ' callingFunction ':' newline];
  fprintf(LOGGER_LOGS{nlog}, str);
  fprintf(LOGGER_LOGS{nlog}, string);  
  fprintf(LOGGER_LOGS{nlog}, '\n\n');
  

  if print_on_stdout
    fprintf(1, str);
    fprintf(1, string);  
    fprintf(1, '\n\n');
  end % if print_on_stdout
  
  
    