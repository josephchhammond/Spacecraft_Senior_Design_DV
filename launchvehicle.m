function [m_launch] = launchvehicle(v_inf)
% Given a v_inf, calculate the Launch Mass for the Falcon Heavy Expendable
% 
% Inputs: 
%   V_inf = Hyperbolic excess speed, m/s
% 
% Outputs:
%   m_launch: Maximum Launch Mass, kg
% 
% Written by Max McDermott for AERO 448 Senior Design 


C3 = (v_inf/1000)^2;
n = 20; 
m = zeros(1,n); 
m(1) = 500; 


for i = 1:n
    
    m(:,i+1) = i*1000; 
    
end

m_vec = m'; 
C3_vec = [12.16552506;11.7046999;10.90871211;10.19803903;9.539;8.944;8.485;7.874;7.416;7;6.403;6.08;5.5677;5.196152423;4.69041576;4.123105626;3.741657387;3.16227766;2.645751311;2;1];

%Polynomial Fit ~ R^2 = 0.99 (Best Fit) 
p = polyfit(C3_vec,m_vec,2); 
p1 = p(1); 
p2 = p(2);
p3 = p(3); 


m_launch = p1*C3^2 + p2*C3 + p3;

end

