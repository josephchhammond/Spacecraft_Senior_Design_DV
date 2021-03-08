function [x1,x2] = quad(a,b,c)
%Quadratic equation solution

const1 = -b./(2.*a);
const2 = sqrt(b.^2 - 4.*a.*c)./(2.*a);

x1 = const1 + const2;
x2 = const1 - const2;
end

