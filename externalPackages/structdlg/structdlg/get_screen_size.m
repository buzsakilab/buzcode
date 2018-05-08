function siz = get_screen_size(units)
%

% AF 10/24/01

if (exist('units','var') ~= 1)
   units = 'pixels';
end

old_units = get(0,'units');
set(0,'units',units);
siz = get(0,'ScreenSize');
set(0,'units',old_units);
