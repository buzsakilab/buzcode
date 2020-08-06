function []=v2structDemo2()
%  V2STRUCTDEMO2 is a usage example demo function for the function V2STRUCT. This demo
%  shows a suggestion for how to use v2struct in managing input to other functions with
%  improved usability.
%
%% Last update:
%    30.06.2011, Adi N.

fprintf(['Running usage demo function for vs2truct function:\n',...
         '\nThis demo shows a suggestion for how to use v2struct in managing input to',...
         ' other functions\nwith improved usability. Please read comments in the .m file for explanations\n',...
         'Press any key between breaks.\n']);
pause()

fprintf('\nOutput from four statements:\n');
pause()
% See the differences between calling foo and fooSmart using v2struct
% Note that we expect the same results in the following statements

% regular call to foo
foo(1,2,3,4);

% pack to structure and call with fooSmart
S = v2struct(1,2,3,4,{'fieldNames','num1','denum1','num2','denum2'});
fooSmart(S);
 
% regular call to fooSmartParamOrder
fooSmartParamOrder(1,2,3,4);

% call with structure to fooSmartParamOrder
fooSmartParamOrder(S);


% Consider the following functions

% foo is a function with regular construction using arguments passed in a defined order
% the result is that you have to call it with 4 arguments, and you must remember their order.
 function isFirstBiggerThanSecond = foo(num1,denum1,num2,denum2)
 isFirstBiggerThanSecond = (num1/denum1>num2/denum2);
 fprintf('\n using foo with a regular call:\n')
 if isFirstBiggerThanSecond
    fprintf('  %d/%d is bigger than %d/%d\n',num1,denum1,num2,denum2);
 else
    fprintf('  %d/%d is not bigger than %d/%d\n',num1,denum1,num2,denum2);
 end

% Now consider the function fooSmart which is modified just slightly to accept input
% packed in a structure. You get the same result but you can use a structure.
    function isFirstBiggerThansSecond = fooSmart(num1,denum1,num2,denum2)
 if isstruct(num1) % unpack if called with a structure
    v2struct(num1);
 end
 isFirstBiggerThansSecond = (num1/denum1>num2/denum2);
 fprintf('\n using fooSmart with a structure:\n')
 if isFirstBiggerThansSecond
    fprintf('  %d/%d is bigger than %d/%d\n',num1,denum1,num2,denum2);
 else
    fprintf('  %d/%d is not bigger than %d/%d\n',num1,denum1,num2,denum2);
 end

% Lastly the function fooSmartParamOrder represent a situation in which the order of the
% parameters is changed in relation to foo. Here calling regularly with the same variables
% as before leads to an erroneous result. However when using the structure as input the
% order is not important and you get the correct result.
     function isFirstBiggerThansSecond = fooSmartParamOrder(num1,num2,denum1,denum2)
 if isstruct(num1) % unpack if called with a structure
    v2struct(num1);
    fprintf('\n using fooSmartParamOrder with structure:\n')
 else
     fprintf('\n using fooSmartParamOrder with same variables:\n')
 end
 isFirstBiggerThansSecond = (num1/denum1>num2/denum2);
 if isFirstBiggerThansSecond
    fprintf('  %d/%d is bigger than %d/%d\n',num1,denum1,num2,denum2);
 else
    fprintf('  %d/%d is not bigger than %d/%d\n',num1,denum1,num2,denum2);
 end