function [rmse] = rmse(value_idea,est_value)
% This function is used to calculate the rmse value
% value_idea: Ideal value
% est_v: filter estimated value

rmse = sqrt(mean((value_idea - est_value).^2));

end

