function s=pretty_print_scientific(x,n)
% pretty-print a value in scientific notation
% usage: s=pp(x,n)
% where: x a floating-point value
%        n is the number of decimal places desired
%        s is the string representation of x with n decimal places, and
%          exponent k
if(x ~= 0)
  exponent=floor(log10(abs(x))); %to accomodate for negative values
  mantissa=x/(10^exponent);
  s=sprintf('$%*.*f \\times 10^{%d}$',n+3,n,mantissa,exponent);
else
  s = sprintf('0');
end
% returns something like '$1.42 \times 10^{-1}$'
end