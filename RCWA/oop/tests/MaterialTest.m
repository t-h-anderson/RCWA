classdef MaterialTest < matlab.unittest.TestCase
    methods (Test)
        function constantMaterial(tc)
            m = rcwa.materials.Constant(3 + 0.5i);
            tc.verifyEqual(m.permittivity([400 500 600]), ...
                (3+0.5i) * ones(1,3), 'AbsTol', 1e-12);
        end

        function glassWithinTableRange(tc)
            m = rcwa.materials.Glass();
            eps = m.permittivity([500 700 1000]);
            tc.verifySize(eps, [1 3]);
            % Real part of n^2 for glass should be ~3.2 (n ~ 1.79)
            tc.verifyGreaterThan(real(eps(1)), 3.0);
            tc.verifyLessThan(real(eps(1)), 3.4);
        end

        function glassNoAbsorptionMatchesLegacy(tc)
            m = rcwa.materials.Glass('absorption', 0);
            lam = [500 700];
            n   = glass('nmlambda', lam);
            tc.verifyEqual(m.permittivity(lam), n.^2, 'RelTol', 1e-12);
        end

        function silverIsComplex(tc)
            m = rcwa.materials.Silver();
            eps = m.permittivity(500);
            tc.verifyTrue(~isreal(eps));
        end

        function aSiHReturnsComplex(tc)
            m = rcwa.materials.ASiH('eVEg', [1.6; 1.7]);
            eps = m.permittivity(600);
            tc.verifySize(eps, [2 1]);
            tc.verifyTrue(any(imag(eps) ~= 0));
        end

        function refractiveIndexConsistent(tc)
            m = rcwa.materials.Constant(4 + 0.01i);
            n = m.refractiveIndex(400);
            tc.verifyEqual(n^2, 4+0.01i, 'AbsTol', 1e-12);
        end
    end
end
