function [Pp,Pv,Kp,Kv,p,v] = KF_standard(n,l_obs,A,x_ini,P_ini,phi,R,Q)
% standard KF

% Input:
% n: the length of the simulation time
% l_obs: observations
% A: observation matrix
% x_ini: initial state vector
% P_ini: initial error matrix
% phi: state transition matrix
% R: observation noise matrix
% Q：system noise matrix

% Output: 
% Pp: P-Matrix(1,1),position-term 
% Pv: P-Matrix(2,2),velocity-term
% Kp: K-Vector(1,1),position-term
% Kv: K-Vector(1,1),velocity-term
% p:position; v:velocity 

%initialization
x_plus_standard = x_ini;
P_plus_standard = P_ini;
for i = 1:n-1    
      Pminus = phi * P_plus_standard * phi' + Q;
      K_standard= inv(A * Pminus * A' + R);
      K_standard = (Pminus * A') * K_standard;
      xhatminus = phi * x_plus_standard;
      x_plus_standard = xhatminus + K_standard * (l_obs(i) - A * xhatminus);
      P_plus_standard = Pminus - K_standard * A * Pminus;
         
      sigma_p(i) = sqrt(P_plus_standard(1,1));                             % standard deviation: position
      sigma_v(i) = sqrt(P_plus_standard(2,2));                             % standard deviation: position    
      Kp(i) = K_standard(1,1);                                             % kalman gain: position    
      Kv(i) = K_standard(2,1);                                             % kalman gain: velocity       
      p(i)=x_plus_standard(1);                                             % filter estimated value: position
      v(i)=x_plus_standard(2);                                             % filter estimated value: velocity
end
Pp = sigma_p';                                                             % P matrix: position 
Pv = sigma_v';                                                             % P matrix: velocity
end

