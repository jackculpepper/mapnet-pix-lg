function [stop,user_data] = cb(x,iter,state,user_data,opts)

switch state
    case 'init'
        user_data.f = zeros(1,opts.maxits);
        user_data.x = zeros(length(x),opts.maxits);
        user_data.its = 0;

        fprintf('\rIteration %4d\tf: %-8.6g', iter.it, iter.f);
    case 'iter'
        user_data.f(iter.it) = iter.f;
        user_data.x(:,iter.it) = x;
        user_data.its = user_data.its + 1;

        fprintf('\rIteration %4d\tf: %-8.6g', iter.it, iter.f);
    case 'done'
        user_data.f = user_data.f(1:user_data.its);
        user_data.x = user_data.x(:,1:user_data.its);

        %fprintf('\r');
end

stop = 0;
