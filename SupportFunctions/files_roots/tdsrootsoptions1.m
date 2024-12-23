function options=tdsrootsoptions1

% the fields concern maximal size of the eigenvalue problem;
% maximal number of newton iterations; the roots accuracy after newton correction; 
% note: if the user doesn't specify the values of options,
% the codes will use the defaulted choice for them, i.e.
% options.max_size_eigenvalue_problem=1000;
% options.newton_max_iterations=20;
% options.root_accuracy=1e-8;
options.max_size_eigenvalue_problem=1000;
options.newton_max_iterations=20;
options.root_accuracy=1e-8;

return;