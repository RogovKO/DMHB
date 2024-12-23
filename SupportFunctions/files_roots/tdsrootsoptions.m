function options=tdsrootsoptions
% options=tdsrootsoptions
%
% the fields concern the minimal real part of the characteristic roots; maximal size of the eigenvalue problem;
% maximal number of newton iterations; the roots accuracy after newton correction; basic delay for commensurate delays; 
%number of grids point p in interval [0, 2*pi) for discretizing \Psi; 
%the value of q which is used for determining new basic delay to
%approximate delays as commensurate delays if m>3 or tau_m/options.commensurate_basic_delay > options.new_basic_delay_q
% note: if the user doesn't specify the values of options,
% the codes will use the defaulted choice for them, i.e.
% options.minimal_real_part=[];
% options.max_size_eigenvalue_problem=1000;
% options.newton_max_iterations=20;
% options.root_accuracy=1e-8;
% options.commensurate_basic_delay=[];
% options.number_grid_points_p=20;
% options.new_basic_delay_q=100;

options.minimal_real_part=[];
options.max_size_eigenvalue_problem=1000;
options.newton_max_iterations=20;
options.root_accuracy=1e-8;
options.commensurate_basic_delay=[];
options.number_grid_points_p=20;
options.new_basic_delay_q=100;
return;