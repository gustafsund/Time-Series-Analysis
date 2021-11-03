function y = myFilter(C,A,x,omitted)

y = filter(C,A,x);
if nargin <4
    omitted = length(A);
end

y = y(omitted:end);

end

