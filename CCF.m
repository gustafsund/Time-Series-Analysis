function [] = CCF(x,y,M)

if nargin <3
    M = 40;
end
stem(-M:M,crosscorr(x,y,M));
title('Cross correlation function')
xlabel('Lags')
hold on
plot(-M:M,2/sqrt(length(x))*ones(1,2*M+1),'-')
plot(-M:M,-2/sqrt(length(x))*ones(1,2*M+1),'-')
hold off

end

