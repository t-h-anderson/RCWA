classdef UnitsTest < matlab.unittest.TestCase
    methods (Test)
        function metresAndNanometres(tc)
            tc.verifyEqual(rcwa.Units.m_from_nm(1000), 1e-6, 'AbsTol', 1e-18);
            tc.verifyEqual(rcwa.Units.nm_from_m(1e-9),   1,   'AbsTol', 1e-12);
        end

        function roundTrip(tc)
            for x = [0.1, 1, 1.234e3, 500]
                tc.verifyEqual(rcwa.Units.nm_from_m(rcwa.Units.m_from_nm(x)), ...
                               x, 'RelTol', 1e-12);
                tc.verifyEqual(rcwa.Units.deg_from_rad(rcwa.Units.rad_from_deg(x)), ...
                               x, 'RelTol', 1e-12);
            end
        end

        function inverseNanometres(tc)
            tc.verifyEqual(rcwa.Units.im_from_inm(1), 1e9, 'AbsTol', 1e-3);
        end

        function degreesAndRadians(tc)
            tc.verifyEqual(rcwa.Units.rad_from_deg(180), pi, 'RelTol', 1e-12);
            tc.verifyEqual(rcwa.Units.deg_from_rad(pi),  180, 'RelTol', 1e-12);
        end

        function photonEnergy(tc)
            % 1.0 eV photon has wavelength 1239.842 nm by definition
            tc.verifyEqual(rcwa.Units.eV_from_nm(1239.842), 1.0, 'RelTol', 1e-9);
            tc.verifyEqual(rcwa.Units.nm_from_eV(1.0), 1239.842, 'RelTol', 1e-9);
        end

        function vectorised(tc)
            v = [400, 500, 600];
            tc.verifyEqual(rcwa.Units.eV_from_nm(v), 1239.842 ./ v, ...
                'RelTol', 1e-12);
        end
    end
end
