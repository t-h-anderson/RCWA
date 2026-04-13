classdef tFourier < matlab.unittest.TestCase
% Tests for Fourier expansion of the grating permittivity.
% Covers epsn.m (analytic) and gepszfft.m (FFT-based) for a rectangular
% grating and checks that both methods agree.

    properties
        loc       % default loc with grating parameters set
        Nt        % number of Fourier modes
        nmLm      % mirror thickness (nm)
        nmLg      % grating relief height (nm)
    end

    methods (TestClassSetup)
        function addSrcToPath(~)
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));
        end
    end

    methods (TestMethodSetup)
        function buildLoc(testCase)
            testCase.Nt   = 2;
            testCase.nmLm = 50;
            testCase.nmLg = 100;

            loc = DefaultLoc;
            loc.Nt   = testCase.Nt;
            loc.nmLx = 400;   % grating period (nm)
            loc.nmLg = testCase.nmLg;
            loc.nmda = testCase.nmLg;  % grating depth = relief (simplest case)
            loc.nmLm = testCase.nmLm;
            loc.zeta = 0.5;   % 50 % duty cycle
            loc.type = 0;     % rectangular grating
            loc.epsm = -10 + 1i;
            loc.epsd = 2.25;
            loc.NFFT = 900;
            testCase.loc = loc;
        end
    end

    methods (Test)

        % ------------------------------------------------------------------
        % Output dimensions
        % ------------------------------------------------------------------

        function testOutputSizeNzBy4Nt1(testCase)
            nmz  = testCase.nmLm + (10:10:90);  % 9 z-positions inside grating
            epsz = epsn(testCase.loc, nmz);
            testCase.verifySize(epsz, [9, 4*testCase.Nt + 1]);
        end

        % ------------------------------------------------------------------
        % DC mode (n=0) equals the average permittivity
        % ------------------------------------------------------------------

        function testGratingDCModeEqualsAveragePermittivity(testCase)
            % For a rectangular grating the n=0 Fourier coefficient is the
            % fill-factor weighted average: epsd*(1-zeta) + epsm*zeta
            loc  = testCase.loc;
            nmz  = loc.nmLm + loc.nmLg/2;  % midpoint of grating region
            epsz = epsn(loc, nmz);

            expected_dc = loc.epsd*(1 - loc.zeta) + loc.epsm*loc.zeta;
            testCase.verifyEqual(epsz(1, 2*testCase.Nt + 1), expected_dc, ...
                'AbsTol', 1e-10);
        end

        function testDutyCycleOneIsAllMetal(testCase)
            % zeta = 1: full metal fill → DC mode = epsm, all AC modes = 0
            loc      = testCase.loc;
            loc.zeta = 1.0;
            nmz      = loc.nmLm + loc.nmLg/2;
            epsz     = epsn(loc, nmz);

            Nt = testCase.Nt;
            testCase.verifyEqual(epsz(1, 2*Nt + 1), loc.epsm, 'AbsTol', 1e-10);
            ac = [epsz(1, 1:2*Nt), epsz(1, 2*Nt+2:end)];
            testCase.verifyEqual(ac, zeros(1, 4*Nt), 'AbsTol', 1e-10);
        end

        function testDutyCycleZeroIsAllDielectric(testCase)
            % zeta = 0: full dielectric fill → DC mode = epsd, all AC = 0
            loc      = testCase.loc;
            loc.zeta = 0.0;
            nmz      = loc.nmLm + loc.nmLg/2;
            epsz     = epsn(loc, nmz);

            Nt = testCase.Nt;
            testCase.verifyEqual(epsz(1, 2*Nt + 1), loc.epsd, 'AbsTol', 1e-10);
            ac = [epsz(1, 1:2*Nt), epsz(1, 2*Nt+2:end)];
            testCase.verifyEqual(ac, zeros(1, 4*Nt), 'AbsTol', 1e-10);
        end

        % ------------------------------------------------------------------
        % Mirror region is homogeneous (no grating modulation)
        % ------------------------------------------------------------------

        function testMirrorRegionDCModeIsEpsm(testCase)
            loc  = testCase.loc;
            nmz  = loc.nmLm / 2;   % strictly inside mirror
            epsz = epsn(loc, nmz);
            testCase.verifyEqual(epsz(1, 2*testCase.Nt + 1), loc.epsm, ...
                'AbsTol', 1e-10);
        end

        function testMirrorRegionACModesAreZero(testCase)
            % In a homogeneous mirror all grating harmonics must vanish
            loc  = testCase.loc;
            Nt   = testCase.Nt;
            nmz  = loc.nmLm / 2;
            epsz = epsn(loc, nmz);

            ac = [epsz(1, 1:2*Nt), epsz(1, 2*Nt+2:end)];
            testCase.verifyEqual(ac, zeros(1, 4*Nt), 'AbsTol', 1e-10);
        end

        % ------------------------------------------------------------------
        % Analytic (epsn) and FFT (gepszfft) methods agree for rect. grating
        % ------------------------------------------------------------------

        function testAnalyticAndFFTAgreeRealPart(testCase)
            loc      = testCase.loc;
            loc.NFFT = 2000;   % fine FFT grid for tight tolerance
            nmz      = loc.nmLm + (5:10:95);  % several z-positions in grating

            epsz_a = epsn(loc, nmz);
            epsz_f = gepszfft(nmz, loc);

            testCase.verifyEqual(real(epsz_a), real(epsz_f), 'AbsTol', 0.05);
        end

        function testAnalyticAndFFTAgreeImagPart(testCase)
            loc      = testCase.loc;
            loc.NFFT = 2000;
            nmz      = loc.nmLm + (5:10:95);

            epsz_a = epsn(loc, nmz);
            epsz_f = gepszfft(nmz, loc);

            testCase.verifyEqual(imag(epsz_a), imag(epsz_f), 'AbsTol', 0.05);
        end

    end
end
