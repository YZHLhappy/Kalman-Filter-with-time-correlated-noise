function [u_ar4,l_obs_ar4,R] = Observation_with_AR4Noise_withoutQ(n,phi_AR4,v,sigma_w,phi,H)
% This function is used to generate observations with ar2 noise(without system noise)

% Input:
% n: the length of the simulation time
% phi_AR4: time-correlated noise coefficient
% v: velocity 
% sigma_w: standard deviation of white noise
% wn: system noise
% phi: state transition matrix
% H: observation matrix

% Output:
% l_obs: observations with ar4 noise
% R: observation error matrix
% u_ar4: ar1 noise

c = 0;                                                                     % ar4 model constant term
x_minus = [0;v];
y = sigma_w * randn(n,1);                                                  % white noise term
u(1)=0;       
u(2)=0;  
u(3)=0;
u(4)=0;
for t=5:n
   u(t) = c + phi_AR4(1)*u(t-1) + phi_AR4(2)*u(t-2)+ ...                   % ar4 noise
       phi_AR4(3)*u(t-3)+ phi_AR4(4)*u(t-4) + y(t); 
end
for i = 1:LE_simu
    x_minus = phi * x_minus;
    l_obs(i) = H * x_minus + u(i);                                         
end
R = var(y);                                                                % observation error matrix
u_ar4 = u';
l_obs_ar4 = l_obs';                                                        % observation with ar4 noise
end