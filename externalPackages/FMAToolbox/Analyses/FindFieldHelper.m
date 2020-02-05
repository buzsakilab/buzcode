function [field] =  FindFieldHelper(map,x,y,threshold,circX,circY)
%[field] =  FindFieldHelper(map,x,y,threshold,circX,circY)

%  Helper function to pass arguments to FindField in order to avoid moving
%  around compiled function

   field = FindField(map,x,y,threshold,circX,circY);

end

