% runTests  Run the full RCWA test suite.
%
% Usage (from the RCWA/ directory):
%   runTests          % runs all suites, prints summary
%   runTests verbose  % verbose output per test
%
% Individual suites can also be run directly, e.g.:
%   runtests('tests/tMaterials')
%   runtests('tests/tFourier')
%   runtests('tests/tRCWA')
%   runtests('tests/tValidation')

addpath(genpath(fullfile(fileparts(mfilename('fullpath')))));

suites = { ...
    'tests/tMaterials',  ...
    'tests/tFourier',    ...
    'tests/tRCWA',       ...
    'tests/tValidation', ...
};

results = cellfun(@(s) runtests(s), suites, 'UniformOutput', false);
results = vertcat(results{:});

table(results)
fprintf('\nPassed: %d / %d\n', sum([results.Passed]), numel(results));
