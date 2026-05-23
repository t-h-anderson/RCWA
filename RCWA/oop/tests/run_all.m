function results = run_all()
% RUN_ALL  Execute the full OOP RCWA test suite.
%
%   results = run_all() returns a TestResult array. Call with no
%   output to just print the summary.
%
% Requires the original RCWA/ directory (and Materials/) on the path
% because several OOP classes wrap legacy tabulated material data.

% Resolve paths
thisDir = fileparts(mfilename('fullpath'));
oopDir  = fileparts(thisDir);
rcwaDir = fileparts(oopDir);

addpath(genpath(rcwaDir));     % original code + materials
addpath(oopDir);               % packages live one above tests/

import matlab.unittest.TestSuite;
suite = TestSuite.fromFolder(thisDir);
results = run(suite);

if nargout == 0
    disp(results);
end
end
