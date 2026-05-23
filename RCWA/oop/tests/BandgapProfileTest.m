classdef BandgapProfileTest < matlab.unittest.TestCase
    methods (Test)
        function legacyReproducesOriginal(tc)
            % With legacy=true (default), match EgProfile.m exactly.
            cfg = rcwa.Config('Eg0', 1.4, 'A', 0.3, 'kappa', 2, ...
                              'phi', 0.1, 'alpha', 2, ...
                              'nmLp', 10, 'nmLi', 200, 'nmLn', 10);
            bp = rcwa.BandgapProfile.fromConfig(cfg);

            nmz = linspace(0, 220, 51);
            Eg_oop    = bp.evaluate(nmz);
            Eg_legacy = EgProfile(nmz, cfg.toStruct());
            tc.verifyEqual(Eg_oop, Eg_legacy, 'AbsTol', 1e-12);
        end

        function pAndNLayersAreWindowGap(tc)
            bp = rcwa.BandgapProfile('Eg0', 1.6, 'A', 0, 'nmLp', 10, ...
                                     'nmLi', 50, 'nmLn', 10);
            % Sample squarely inside p-, i-, n- regions
            Eg = bp.evaluate([5, 40, 65]);
            tc.verifyEqual(Eg(1), 1.95, 'AbsTol', 1e-12);
            tc.verifyEqual(Eg(2), 1.6,  'AbsTol', 1e-12);
            tc.verifyEqual(Eg(3), 1.95, 'AbsTol', 1e-12);
        end

        function flatProfileWhenAmplitudeZero(tc)
            bp = rcwa.BandgapProfile('Eg0', 1.7, 'A', 0, ...
                                     'nmLp', 0, 'nmLi', 100, 'nmLn', 0, ...
                                     'legacy', false);
            % Sample strictly inside the i-layer (0,100). The boundaries
            % land on the layer-categorisation sign() seams and fall into
            % no layer (faithful to the original sign-trick semantics).
            nmz = linspace(1, 99, 21);
            tc.verifyTrue(all(abs(bp.evaluate(nmz) - 1.7) < 1e-12));
        end

        function bugDocumentedRegression(tc)
            % BUG-001: legacy offsets the i-layer sin by nmLn rather
            % than nmLp. Demonstrate that when nmLp != nmLn the legacy
            % and corrected profiles diverge.
            cfg = rcwa.Config('Eg0', 1.5, 'A', 0.3, 'kappa', 3, ...
                              'phi', 0, 'alpha', 1, ...
                              'nmLp', 20, 'nmLi', 200, 'nmLn', 5);
            bp_legacy = rcwa.BandgapProfile.fromConfig(cfg);
            bp_fixed  = bp_legacy;
            bp_fixed.legacy = false;
            nmz = linspace(20, 220, 41);
            diff = bp_legacy.evaluate(nmz) - bp_fixed.evaluate(nmz);
            tc.verifyGreaterThan(max(abs(diff)), 1e-3, ...
                'Expected legacy vs corrected to differ when nmLp != nmLn');
        end
    end
end
