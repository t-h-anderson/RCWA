classdef tRCWA < matlab.unittest.TestCase
% Tests for the RCWA engine (RCWA.m) and specular coefficient extraction
% (specular.m).
%
% All tests use a single lossless glass slab (n=1.5, eps=2.25) in air at
% normal incidence. The grating period is chosen so that all diffraction
% orders except the zeroth are evanescent, making the system equivalent to
% a standard Fabry-Pérot etalon with analytic solutions.
%
% Geometry:
%   incident medium (air, nsa=1)  |  glass slab (eps=2.25)  |  air (nsa=1)
%
% Parameters:
%   lambda  = 600 nm  → k0 = 2π/600 nm⁻¹
%   nmLx    = 300 nm  → first-order kx = 2π/300 > k0, so evanescent in air
%   Nt      = 2       → 5 Fourier modes (-2,-1,0,1,2)

    properties
        loc       % pre-configured DefaultLoc
        Nt
        inmk0
        radtheta
        epsf_unit % epsf for a single homogeneous glass layer
    end

    methods (TestClassSetup)
        function addSrcToPath(~)
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));
        end
    end

    methods (TestMethodSetup)
        function buildInputs(testCase)
            testCase.Nt       = 2;
            testCase.inmk0    = 2*pi/600;  % k0 at 600 nm (nm⁻¹)
            testCase.radtheta = 0;         % normal incidence

            loc = DefaultLoc;
            loc.Nt        = testCase.Nt;
            loc.nmLx      = 300;   % nm; makes n=±1 orders evanescent at 600 nm
            loc.pol       = 1;     % p-polarisation
            loc.nsa       = 1.0;   % exact real air (avoids tiny energy leak from default 1+1e-6i)
            loc.parforArg = 0;
            testCase.loc = loc;

            % Fourier coefficients for a homogeneous glass layer (eps=2.25):
            % only the zero-order term (column 2*Nt+1) is non-zero.
            Nt = testCase.Nt;
            ef = zeros(1, 4*Nt + 1);
            ef(1, 2*Nt + 1) = 2.25;
            testCase.epsf_unit = ef;
        end
    end

    % ----------------------------------------------------------------------
    % Helper
    % ----------------------------------------------------------------------
    methods (Access = private)
        function [rp0, rp1, tp0, tp1] = runPlanar(testCase, d_nm)
            % Run RCWA + specular for a single glass layer of thickness d_nm.
            [~, Tn, R, ~] = RCWA(testCase.inmk0, testCase.radtheta, ...
                testCase.epsf_unit, d_nm, testCase.loc);
            [rp0, rp1, tp0, tp1] = specular(R, Tn{end}, ...
                testCase.inmk0, testCase.radtheta, testCase.loc);
        end
    end

    methods (Test)

        % ------------------------------------------------------------------
        % Output dimensions
        % ------------------------------------------------------------------

        function testOutputDimensionsR(testCase)
            % R must be a column vector of length 2*Nt+1
            [~, ~, R, ~] = RCWA(testCase.inmk0, testCase.radtheta, ...
                testCase.epsf_unit, 100, testCase.loc);
            testCase.verifySize(R, [2*testCase.Nt + 1, 1]);
        end

        function testOutputDimensionsTn(testCase)
            % Tn is a cell array of length Ns+1 = 2 for a single-layer stack
            [~, Tn, ~, ~] = RCWA(testCase.inmk0, testCase.radtheta, ...
                testCase.epsf_unit, 100, testCase.loc);
            testCase.verifySize(Tn, [2, 1]);
            testCase.verifySize(Tn{end}, [2*testCase.Nt + 1, 1]);
        end

        function testOutputDimensionsE(testCase)
            % E has shape (3, Ns, 2*Nt+1) = (3, 1, 5) for a single layer
            [~, ~, ~, E] = RCWA(testCase.inmk0, testCase.radtheta, ...
                testCase.epsf_unit, 100, testCase.loc);
            testCase.verifySize(E, [3, 1, 2*testCase.Nt + 1]);
        end

        % ------------------------------------------------------------------
        % Energy conservation: R + T = 1 for a lossless slab
        %
        % With nsa=1 exactly and real eps=2.25, and all higher diffraction
        % orders evanescent, the specular power must satisfy rp0+tp0 = 1.
        % ------------------------------------------------------------------

        function testEnergyConservationLosslessSlab(testCase)
            [rp0, ~, tp0, ~] = testCase.runPlanar(100);
            testCase.verifyEqual(rp0 + tp0, 1.0, 'AbsTol', 1e-6, ...
                'R + T should equal 1 for a lossless slab in air');
        end

        function testEacvanascentOrdersCarryNoPower(testCase)
            % rp1 and tp1 come from evanescent orders (Re(kz)=0), so both
            % should be negligibly small
            [~, rp1, ~, tp1] = testCase.runPlanar(100);
            testCase.verifyLessThan(abs(rp1), 1e-10);
            testCase.verifyLessThan(abs(tp1), 1e-10);
        end

        % ------------------------------------------------------------------
        % Half-wave condition: R = 0
        %
        % Glass slab, d = 200 nm, lambda = 600 nm:
        %   phi = k0 * n * d = (2π/600) * 1.5 * 200 = π  (half-wave)
        % Fabry-Pérot: r = (r01 + r12*exp(2iπ)) / (1 + r01*r12*exp(2iπ))
        %                 = (−0.2 + 0.2) / (1 − 0.04) = 0
        % ------------------------------------------------------------------

        function testHalfWaveZeroReflectance(testCase)
            [rp0, ~, ~, ~] = testCase.runPlanar(200);
            testCase.verifyEqual(rp0, 0, 'AbsTol', 1e-10, ...
                'Half-wave slab should give zero specular reflectance');
        end

        function testHalfWaveFullTransmittance(testCase)
            [~, ~, tp0, ~] = testCase.runPlanar(200);
            testCase.verifyEqual(tp0, 1.0, 'AbsTol', 1e-10, ...
                'Half-wave slab should give unity specular transmittance');
        end

        % ------------------------------------------------------------------
        % Absorbing slab: R + T < 1
        % ------------------------------------------------------------------

        function testAbsorbingSlabDissipatesEnergy(testCase)
            Nt = testCase.Nt;
            epsf_abs = zeros(1, 4*Nt + 1);
            epsf_abs(1, 2*Nt + 1) = 2.25 + 0.5i;  % lossy glass

            [~, Tn, R, ~] = RCWA(testCase.inmk0, testCase.radtheta, ...
                epsf_abs, 200, testCase.loc);
            [rp0, ~, tp0, ~] = specular(R, Tn{end}, ...
                testCase.inmk0, testCase.radtheta, testCase.loc);

            testCase.verifyLessThan(rp0 + tp0, 1.0, ...
                'Lossy slab must dissipate energy: R + T < 1');
        end

    end
end
