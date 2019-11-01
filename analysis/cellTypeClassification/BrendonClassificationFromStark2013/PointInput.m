function [x,y,button, key] = PointInput(varargin)
%function [x,y,button, key] = PointInput(n,figh, axish)
% Input:
% n = number of clicks
% figh - which figure to analyze the clicks in
% axish - in which axes .. try to avoid this complications, use defaults
% Output:
% x,y - coordinates of the axes where the click happened,
% button - number of button presseed (1 -left,2 -middle, 3 -right
% if key not button was pressed then button=0, and key returns the key that
% was pressed
% also if it is not a character but a system key (like arrow)
% you can get double(key) and use my table from ~antsiro/matlab/keycodes
% file . enjoy :))
% Anton
[n, figh, axish] = DefaultArgs(varargin,{1,gcf,gca});
x=[];y=[];key={};
if n>0
    for i=1:n
        res = waitforbuttonpress;
        if res==0
            whatbutton = get(figh,'SelectionType');
            mousecoord = get(axish,'CurrentPoint');
            x(i)=mousecoord(1,1);
            y(i) = mousecoord(1,2);

            switch whatbutton
                case 'normal'
                    button(i)=1;
                case 'extend'
                    button(i)=2;
                case 'alt'
                    button(i)=3;
            end
            key{end+1}= '';
        else
            x=NaN; y=NaN; button=0;
            key{end+1} = get(gcf,'CurrentCharacter');
        end
    end
    if n==1
        key=key{1};
    end
else
    whatbutton = get(figh,'SelectionType');
    mousecoord = get(axish,'CurrentPoint');
    x= mousecoord(1,1);
    y = mousecoord(1,2);

    switch whatbutton
        case 'normal'
            button=1;
        case 'extend'
            button=2;
        case 'alt'
            button=3;
        otherwise 
   %         fprintf('you pressed button %s\n',whatbutton);
            button = 0;
    end
end

