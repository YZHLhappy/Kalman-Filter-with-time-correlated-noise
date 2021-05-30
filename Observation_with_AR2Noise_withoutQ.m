function [u_ar2,l_obs_ar2,R] = Observation_with_AR2Noise_withoutQ(n,phi_AR2,v,sigma_w)
% This function is used to generate observations with ar2 noise(without system noise)

% Input:
% n: the length of the simulation time
% phi_AR2: time-correlated noise coefficient
% v: velocity 
% sigma_w: standard deviation of white noise

% Output:
% l_obs: observations with ar2 noise
% R: observation noise matrix
% u_ar2: ar1 noise
c = 0;                                                                     % ar2 model constant term
u(1)=0;       
u(2)=0;            
y = sigma_w * randn(n,1);                                                  % white noise term
for t=3:n
   u(t) = c + phi_AR2(1)*u(t-1) + phi_AR2(2)*u(t-2) + y(t);                % ar2 noise
end
l_obs_ar2 = v + u';                                                        % observation with ar2 noise
R = var(y);                                                                % observation matrix
u_ar2 = u';
end