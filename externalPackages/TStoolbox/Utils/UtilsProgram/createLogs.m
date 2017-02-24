function createLogs(logs)
  
  global LOGGER_LOGS
  
  
  if ~iscell(logs)
    lg{1} = logs;
    logs = lg;
  end
  
  for i = 1:length(logs)
    lf{i} = fopen(logs{i}, 'w');
  end
  
  LOGGER_LOGS = lf;
  