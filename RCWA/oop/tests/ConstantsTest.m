classdef ConstantsTest < matlab.unittest.TestCase
    methods (Test)
        function maxwellRelation(tc)
            % c = 1/sqrt(mu0 * eps0), within ~0.1% given the rounded c=3e8.
            c_derived = 1 / sqrt(rcwa.Constants.mu0 * rcwa.Constants.eps0);
            tc.verifyEqual(c_derived, rcwa.Constants.c, 'RelTol', 1e-3);
        end

        function impedanceOfFreeSpace(tc)
            % eta0 = sqrt(mu0/eps0)
            eta = sqrt(rcwa.Constants.mu0 / rcwa.Constants.eps0);
            tc.verifyEqual(eta, rcwa.Constants.eta0, 'RelTol', 1e-3);
        end

        function constantsArePositive(tc)
            tc.verifyGreaterThan(rcwa.Constants.q,    0);
            tc.verifyGreaterThan(rcwa.Constants.h,    0);
            tc.verifyGreaterThan(rcwa.Constants.c,    0);
            tc.verifyGreaterThan(rcwa.Constants.eta0, 0);
        end
    end
end
