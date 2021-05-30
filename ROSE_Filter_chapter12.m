function [R_x_11,R_x_22,Q_x_11,Q_x_22,Q_end,R_end,x,P_11,P_22] =...
    ROSE_Filter_chapter12(x_ini,mess,Ts,phi,C)
% This function is the ROSE-Filter algorithm
% from "Kalman Filter",by Reiner Marchthaler·Sebastian Dingler,chapter 12
% Input：
%      x_ini: Initial state,[position;velocity]
%      mess:Measurement matrix [2x360],mess = [position;velocity]
%      Ts: the sampling period
%      phi: state transition matrix
%      C: observation matrix

% Output:
% R: measurement noise matrix
% Q: process noise matrix
% R_x_11: R-matrix(1,1),position term
% R_x_22: R-matrix(2,2),velocity term
% Q_x_11: Q-matrix(1,1),position term
% Q_x_22: Q-matrix(2,2),velocity term
% Q_end: the last Q in the iteration
% R_end: The last R in the iteration
% x: filter estimated value
% P_11: P matrix(1,1),position term
% P_22: P matrix(2,2),velocity term

x_x_plus = x_ini;
x_x_minus = x_ini;
y_x_mess = mess;                                                           % [position;velocity];
%ROSE Parameters
Q_ini = 1;
R_ini = 1;
c_x = [1,0];
lamba = Ts*sqrt(Q_ini/R_ini);                                              % constant term
K_x = (0.125/Ts).*[Ts * (-lamba^2-8*lamba+(lamba+4)*sqrt(lamba^2+8*lamba));% kalman gain for R
                   2*(lamba^2+4*lamba-lamba*sqrt(lamba^2+8*lamba))];
H_x = (eye(length(phi))-K_x * c_x) * phi;
x_x1 = x_ini;
x_x2 = [mess(2,1);0];
Gamma = 2.5;                                                               % Gain factor measurement noise                                                        
Alpha_R = 0.05;                                                            % Kalman-gain variance measurement noise
Alpha_M = 0.1;                                                             % Kalman-gain variance M
R =[1,0;0,1e-1];
P_plus = [1e-1,0;0,1e-2];
M = C * P_plus * C';
Q_max = 3e-1;                                                              % min.Q-Wert
Q_min = 3e-8;                                                              % max. Q-Wert
for i = 1:360
   % Bestimmung R durch ROSE-Filter
   y1_x = y_x_mess(1,i);
   y2_x = y_x_mess(2,i);
   x_x1 = H_x * x_x1 + K_x * y1_x;
   x_x2 = H_x * x_x2 + K_x * y2_x;
   R = Gamma * Alpha_R *[x_x1(1)-y1_x;x_x2(1)-y2_x]*...
       [x_x1(1)-y1_x;x_x2(1)-y2_x]'+(1-Alpha_R)*R;
   R_x_11(i) = R(1,1);
   R_x_22(i) = R(2,2);
   
   % Bestimmung Q durch ROSE-Filter
   dy_x = y_x_mess(:,i)-C*x_x_plus;
   M = Alpha_M .* dy_x * dy_x'+(1-Alpha_M).* M;
   Q = C'*(M - R)*C -phi * P_plus * phi';
   
   if Q(1,1) < Q_min
       Q(1,1) = Q_min;
   end
   if Q(2,2) < Q_min
       Q(2,2) = Q_min;
   end
   if Q(1,1) > Q_max
       Q(1,1)=Q_max;
   end
   if Q(2,2) > Q_max
       Q(2,2)=Q_max;
   end
      Q_x_11(i) = Q(1,1);
      Q_x_22(i) = Q(2,2);
      
   % Kalman Filter part
   p_minus = phi * P_plus *phi +Q;
   K = p_minus*C' * pinv(C*p_minus*C'+R);
   x_x_plus = x_x_minus+K*dy_x;
   P_plus = (eye(length(p_minus))-K*C)*p_minus;
   x_x_minus = phi * x_x_plus;
   
   Q_end = Q;
   R_end = R;
   x(:,i) = x_x_plus;
   P_11(i) = P_plus(1,1);
   P_22(i) = P_plus(2,2);
end
end

