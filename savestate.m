
vars = who;
vars_str = sprintf('%s ', vars{:});

%% don't save Xr or cache
vars_str = strrep(vars_str, ' Xr ', ' ');
vars_str = strrep(vars_str, ' Zr ', ' ');
vars_str = strrep(vars_str, ' cache ', ' ');

eval(sprintf('save state/%s/matlab_up=%06d.mat %s', paramstr, update, vars_str)); 

