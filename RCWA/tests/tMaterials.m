classdef tMaterials < matlab.unittest.TestCase
% Unit tests for material permittivity functions:
% glass, Silver, AZO, aSiHGC

    methods (TestClassSetup)
        function addSrcToPath(~)
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));
        end
    end

    methods (Test)

        % ------------------------------------------------------------------
        % glass
        % ------------------------------------------------------------------

        function testGlassKnownValue550nm(testCase)
            % Hardcoded table value at exactly 550 nm
            val = glass('nmlambda', 550);
            testCase.verifyEqual(val, 1.79115, 'AbsTol', 1e-5);
        end

        function testGlassKnownValue800nm(testCase)
            % Hardcoded table value at exactly 800 nm
            val = glass('nmlambda', 800);
            testCase.verifyEqual(val, 1.76462, 'AbsTol', 1e-5);
        end

        function testGlassReturnsRealValue(testCase)
            % glass models only the real refractive index (no absorption term)
            val = glass('nmlambda', 700);
            testCase.verifyTrue(isreal(val));
        end

        function testGlassNormalDispersion(testCase)
            % Refractive index decreases monotonically with wavelength
            lam  = [500, 700, 1000];
            vals = arrayfun(@(l) glass('nmlambda', l), lam);
            testCase.verifyTrue(all(diff(vals) < 0));
        end

        function testGlassVectorInput(testCase)
            lam = [500, 700, 900];
            val = glass('nmlambda', lam);
            testCase.verifySize(val, size(lam));
        end

        % ------------------------------------------------------------------
        % Silver
        % ------------------------------------------------------------------

        function testSilverScalarOutput(testCase)
            val = Silver('nmlambda', 600);
            testCase.verifySize(val, [1, 1]);
        end

        function testSilverIsComplex(testCase)
            % Refractive index of silver has both n and k components
            val = Silver('nmlambda', 600);
            testCase.verifyFalse(isreal(val));
        end

        function testSilverPositiveExtinctionCoefficient(testCase)
            % k > 0: silver absorbs at visible wavelengths
            val = Silver('nmlambda', 400);
            testCase.verifyGreaterThan(imag(val), 0);
        end

        function testSilverPositiveRealIndex(testCase)
            val = Silver('nmlambda', 400);
            testCase.verifyGreaterThan(real(val), 0);
        end

        function testSilverVectorInput(testCase)
            lam = [400, 600, 800];
            val = Silver('nmlambda', lam);
            testCase.verifySize(val, [1, 3]);
        end

        % ------------------------------------------------------------------
        % AZO (aluminium-doped zinc oxide)
        % ------------------------------------------------------------------

        function testAZOScalarOutput(testCase)
            val = AZO(600);
            testCase.verifySize(val, [1, 1]);
        end

        function testAZOPositiveRealIndex(testCase)
            val = AZO(600);
            testCase.verifyGreaterThan(real(val), 0);
        end

        function testAZONonNegativeExtinctionCoefficient(testCase)
            % k >= 0 (passive material)
            val = AZO(600);
            testCase.verifyGreaterThanOrEqual(imag(val), 0);
        end

        function testAZOVectorInput(testCase)
            lam = [400, 600, 800];
            val = AZO(lam);
            testCase.verifySize(val, [1, 3]);
        end

        % ------------------------------------------------------------------
        % aSiHGC (Cody-Lorentz model for amorphous silicon)
        % ------------------------------------------------------------------

        function testaSiHGCSingleEgSingleLambdaShape(testCase)
            out = aSiHGC(600, 1.6);
            testCase.verifySize(out, [1, 1]);
        end

        function testaSiHGCMultiEgSingleLambdaShape(testCase)
            Eg  = [1.5; 1.6; 1.7; 1.8];
            out = aSiHGC(600, Eg);
            testCase.verifySize(out, [4, 1]);
        end

        function testaSiHGCMultiLambdaShape(testCase)
            Eg       = [1.5; 1.7];
            nmlambda = [400, 600, 800];
            out = aSiHGC(nmlambda, Eg);
            testCase.verifySize(out, [2, 3]);
        end

        function testaSiHGCNonNegativeImagPart(testCase)
            % eps2 >= 0 everywhere (passive material, no gain)
            Eg       = linspace(1.4, 2.0, 5)';
            nmlambda = [400, 600, 800, 1000];
            out = aSiHGC(nmlambda, Eg);
            testCase.verifyGreaterThanOrEqual(imag(out), zeros(5, 4));
        end

        function testaSiHGCAboveGapIsAbsorbing(testCase)
            % 400 nm = 3.1 eV >> 1.6 eV bandgap → eps2 > 0
            out = aSiHGC(400, 1.6);
            testCase.verifyGreaterThan(imag(out), 0);
        end

        function testaSiHGCBelowGapIsTransparent(testCase)
            % 900 nm = 1.38 eV < 1.6 eV bandgap → eps2 ≈ 0
            out = aSiHGC(900, 1.6);
            testCase.verifyLessThan(imag(out), 0.1);
        end

        function testaSiHGCHigherBandgapLessAbsorption(testCase)
            % At fixed wavelength 500 nm (2.48 eV), a higher bandgap
            % material absorbs less (photon energy is "less above" the gap)
            nmlambda = 500;
            eps_low  = aSiHGC(nmlambda, 1.6);
            eps_high = aSiHGC(nmlambda, 2.2);
            testCase.verifyGreaterThan(imag(eps_low), imag(eps_high));
        end

    end
end
