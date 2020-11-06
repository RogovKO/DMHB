function c = chebpoly(n, x)
%CHEBPOLY Chebyshev polynomial of the first kind.
%
%   CHEBPOLY(N) returns the coefficients of the polynomial of degree N.
%
%   CHEBPOLY(N, X) returns the polynomial of degree N evaluated in X.
%
%
     

      if n==-1
          c=zeros(size(x));
      elseif n == 0                 % use explicit formula when N = 0
         c = ones(size(x));
      elseif n == 1             % use explicit formula when N = 1
         c = x;
      else                      % use recursive formula when N > 1
         a = 1;
         b = x;
         for k = 2 : n
            c = 2 * x .* b - a;
            a = b;
            b = c;
         end
      end

     