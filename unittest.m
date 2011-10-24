
clear

J = 31^2;
J = 25^2;
J = 19^2;
J = 15^2;

Lmin = 20;
Lstp = 5;
Lmax = 100;

N = 25;
N = 9;
N = 16;

Jsz = sqrt(J);
Nsz = sqrt(N);


display_every = 1000;
save_every = 1000;
c_log_length = 1000;
b_log_length = 1000;

datasource = 'movies'; datatype = 'chunks';

lambda = 0.1;
lower_var_thresh = 0.01;

inf_type = 'lbfgsb';
lrn_type = 'gd';

opts_lbfgs_c = lbfgs_options('iprint', -1, 'maxits', 20, ...
                           'factr', 0.01, ...
                           'cb', @cb);


paramstr = sprintf('J=%03d_L=%03d_N=%03d_%s',J,Lmax,N,datestr(now,30));

B = 10000;
Bmini = 10;
eta = 0.005;
eta = 0.1;
eta = 0.03;
eta = 0.001;

c_log = zeros(N, Bmini, c_log_length);
b_log = zeros(2, Bmini, b_log_length);

target_angle = 0.005;
angle_thresh = 0.015;

reinit;


border = 5;
border = 2;
mask = zeros(Jsz,Jsz);
valid_range = border+1:Jsz-border;
mask(valid_range,valid_range) = ones(length(valid_range));
mask = reshape(mask, J, 1);



num_trials = 5000;
num_trials = 2000;
num_trials = 10000;

L = Lmin;

for L = L:Lstp:Lmax ; mapnet ; end

