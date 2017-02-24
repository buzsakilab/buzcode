  lv = length(varargin);
  
  i = 1;
  done = 0;
  while ~done
    cs = ['%' num2str(i)];
    vs = ['varargin{' num2str(i) '}'];
    if isempty(findstr(cmd_string, cs))
      done = 1;
      break;
    end
    cmd_string = substr(cmd_string, cs, vs);
  end