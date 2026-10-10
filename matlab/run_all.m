%RUN_ALL Regenerate everything the RTL flow needs from the MATLAB model.
%   1) ROM .mem files under rtl/   2) model self-checks   3) testbench vectors
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'model'));

gen_luts();
assert(check_model(false), 'Golden model self-check failed');
gen_vectors();
