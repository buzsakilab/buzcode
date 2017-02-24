function closeLogs()
  
  global LOGGER_LOGS

  lf = LOGGER_LOGS;
  
  for i = 1:length(lf)
    fclose(lf{i});
  end
  