function [u_ar1,l_obs,R] = Observation_with_AR1Noise_withoutQ(n,psi,v,sigma_ce,Ts,tau)
% This function is used to generate observations with ar1 noise(without system noise)  
% (hier ar1: Gaussian Markov Model)

% Input:
% n: the length of the simulation time
% psi: time-correlated noise coefficient
% v: velocity 
% tau: time constant
% Ts: period [s]
% sigma_ce: standard deviation of correlated errors

% Output:
% l_obs: observations with ar1 noise
% R: observation noise matrix
% u_ar1: ar1 noise  

sigma_W_GM = sqrt(sigma_ce^2 * (1-psi^2));                                 % standard deviation of white noise
W_GM = sigma_W_GM*randn(LE_simu,1);                                        % white noise
u = 0;                                                                     % initialize u
for j=1:n-1
    u(j+1) = u(j) * exp(-Ts/tau) + W_GM(j);                                % u is time-correlated error
end
l_obs = v+u';                                                              % observation
R = var(W_GM);                                                             % observation noise matrix
u_ar1 = u';                                                                % ar1 noise                                                 
end
