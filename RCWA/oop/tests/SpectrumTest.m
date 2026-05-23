classdef SpectrumTest < matlab.unittest.TestCase
    methods (Test)
        function matchesLegacyAtTablePoints(tc)
            spec = rcwa.Spectrum();
            for nm = [400, 500, 700, 1000]
                tc.verifyEqual(spec.spectralIrradiance(nm), ...
                               W2inmim2AM15G('nmlambda', nm), ...
                               'RelTol', 1e-12);
            end
        end

        function interpolatesBetweenSamples(tc)
            spec = rcwa.Spectrum();
            v = spec.spectralIrradiance([500.25, 500.75]);
            tc.verifyTrue(all(isfinite(v)));
        end
    end
end
