function [check_S] = check_convergence(vector_S,Sk,size)
% Get if the condition is met
% condition_S = zeros(size,1);
% order_S = zeros(size,1);
% for j = 1:size
%     if floor(log10(vector_S(j))) == 0 || vector_S(j) < 10^(-15)
%         order_S(j) = 1;
%     else
%         order_S(j) = floor(log10(vector_S(j)));
%     end
%     condition_S(j) = gt(10^(order_S(j) - 5),abs(vector_S(j) - Sk(j)));
% end
% check_S = not(all(condition_S == 1));
rms = sqrt(sum((vector_S - Sk).^2)/size);
if rms <= 1e-4
    check_S = false;
else
    check_S = true;
end
end

