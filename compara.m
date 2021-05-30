function [a_better] = compara(rmse_a,rmse_b,j)
% This function is used to compare the rmse value, and output the number 
% of times the rmse value of a is smaller

% Input:
% rmse_a: the rmse value of the a
% rmse_b: the rmse value of the b
% the number of rmse value

% Output:
% the number of times the rmse value of a is smaller

diff = rmse_b - rmse_a;
[~,b] = size(diff);
b_good_than_a = 0;
for f = 1:b
   if diff(f) < 0
      b_good_than_a = b_good_than_a +1;
   end
end
 a_better = j - b_good_than_a;   
end

