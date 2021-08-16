function i = issame2(x,y)
% i = issame2(x,y)
%
% Tests whether the values of the variables X and Y are the same.  For 
% arbitrary matrices X and Y, ISSAME(X,Y) returns I=1 if 
% double(X)==double(Y).  It returns I=2 if, additionally,
% ISSTR(X)==ISSTR(Y).  Otherwise, it returns I=0.
%

i=0;
if isequal(double(x),double(y)),
   if ~abs(max(isnan(x(:))-isnan(y(:)))),
       i=1;
       if isstr(x)==isstr(y),
          i=2;
       end
   end
end
