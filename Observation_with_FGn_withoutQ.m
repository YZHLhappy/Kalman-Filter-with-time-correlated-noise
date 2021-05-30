function [u_out,l_obs_FGn,R] = Observation_with_FGn_withoutQ(n_simu,v,A,phi,sigma,H,n,N,force)
% This function is used to generate observations with Fractional Gaussian noise(without system noise) 

% Input:
% n_simu: the length of the simulation time
% v: velocity 
% A: observation matrix
% phi: state transition matrix
% sigma,H,n,N,force are the parameter of ffgn function.

% Output:
% l_obs_FGn: observations with Fractional Gaussian noise
% R: observation noise matrix
% u_out: Fractional Gaussian noise  

x_minus = [0;v];                                                           % initialize state vector
u = ffgn(sigma,H,n,N,force);                                               % Fractional Gaussian noise

for i = 1:n_simu
    x_minus = phi * x_minus;
    l_obs(i) = A * x_minus + u(i);
end
R = var(u);                                                                % observation noise matrix
u_out = u';                                                                
l_obs_FGn = l_obs';                                                        % observation
end

