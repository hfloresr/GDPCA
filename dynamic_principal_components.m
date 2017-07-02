% This is a simple implementation of the dynamic principal component anaylsis
% methods presented in the following paper.
%
% Pe�a, Daniel, and Victor J. Yohai. "Generalized Dynamic Principal Components."
% Journal of the American Statistical Association 111.515 (2016): 1121-1131.
%
% TODOs: 
% - (1.) vectorize some operations since current implementation using
% may explicit for loops which can be extremely inefficient
% - (2.) refactor the code and wrap them in a class for this model. 
% The class may contains methods like model.train(), model.evaluate(), ...
% - (3.) add test cases (currently zero test case)

function dynamic_principal_components
% run the DPC model

% set up
%z = normrnd(0, 1, [10, 100]);
z = load('z.mat')
z = z.z
k1 = 5;
k2 = 3;
k = k1 + k2;
[m, T] = size(z);
%f_init = normrnd(0, 1, [T + k, 1]);
f_init = load('f.mat')
f_init = f_init.f
alpha_init = nan;
beta_init = nan;

% train
f = f_init;
alpha = alpha_init;
beta = beta_init;
train_iterations = 100;
loss_values = zeros(train_iterations, 1);

% Save z and f for testing in R
%save('z.mat', 'z')
%save('f.mat', 'f')

for train_iter = 1:train_iterations
[f, alpha, beta] = run_train_step(z, k, f, alpha, beta);
loss_values(train_iter) = evaluate_op(z, k, f, alpha, beta);
end

% plot result
figure;
plot(loss_values, 'LineWidth', 1.5);
title('train loss (reconstruction MSE)');
xlabel('iteration');

figure;
plot(f, 'LineWidth', 1.5);
title('factor');
xlabel('time');

% Save loss_values
save('loss_values.mat', 'loss_values')

end

% ---------- below are all helper functions -------------------------------
% not optimized for speed yet (use a lot of for loops)
%
function [f_new, alpha_new, beta_new] = run_train_step(z, k, f, alpha, beta)
[m, T] = size(z);
[alpha_new, beta_new] = alpha_beta(z, f, k);
f_star = f_alpha_beta(z, k, alpha_new, beta_new);
f_centered = f_star -  mean(f_star);
f_new = sqrt(T + k) * f_centered / norm(f_centered);
end

function loss = evaluate_op(z, k, f, alpha, beta)
loss = mean_squared_error(z, k, f, alpha, beta);
end

function C = C_alpha(z, alpha, k)
[m, T] = size(z);
C = zeros(m, T + k, k + 1);
for j = 1:m
    for t = 1:size(C, 2)
        for q = 1:size(C, 3)
            if (q >= max(1, t - T + 1)) && (q <= min(k + 1, t))
                C(j, t, q) = z(j, t - q + 1) - alpha(j);
            end
        end
    end
end
end

function D = D_beta(z, beta, k)
[m, T] = size(z);
D = zeros(m, T + k, T + k);
for j = 1:m
    for t = 1:size(D, 2)
        for q = 1:size(D, 3)
            if (q >= max(t - k, 1)) && (q <= min(t + k, T + k))
                for v = max([1, t - k, q - k]):min([t, q, T])
                    D(j, t, q) = D(j, t, q) + beta(j, q - v + 1) * beta(j, t - v + 1);
                end
            end
        end
    end
end

end

function f = f_alpha_beta(z, k, alpha, beta)
[m, T] = size(z);
D_3d = D_beta(z, beta, k);
D = squeeze(sum(D_3d, 1));
C_3d = C_alpha(z, alpha, k);
sum_Cj_betaj = 0;
for j = 1:m
    sum_Cj_betaj = sum_Cj_betaj + squeeze(C_3d(j, :, :)) * beta(j, :)';
end

f = D \ sum_Cj_betaj;
end

function F = F_f(z, f, k)
[m, T] = size(z);
F = zeros(T, k + 2);

for t = 1:size(F, 1)
    F(t, :) = [f(t:(t + k))' 1];
end
end

function [alpha, beta] = alpha_beta(z, f, k)
[m, T] = size(z);
F = F_f(z, f, k);
FtF_inv = inv(F' * F);
FtF_inv_Ft = FtF_inv * F';
tmp = z * FtF_inv_Ft';
alpha = tmp(:, k + 2)';
beta = tmp(:, 1:(k + 1));
end

function err = mean_squared_error(z, k, f, alpha, beta)
[m, T] = size(z);
sum_squared_error = 0;
for j = 1:m
    for t = 1:T
        z_jt_predict = alpha(j);
        for i = 0:k
            z_jt_predict = z_jt_predict + beta(j, i + 1) * f(t + i);
        end
        sum_squared_error = sum_squared_error + (z(j, t) - z_jt_predict)^2;
    end
end
err = sum_squared_error / (T * m);
end