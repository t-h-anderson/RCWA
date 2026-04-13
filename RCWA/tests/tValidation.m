classdef tValidation < matlab.unittest.TestCase
% Tests for the arguments-block input validation added to core functions.
% Each test confirms that an invalid argument is rejected with the expected
% error identifier before any computation begins.

    properties
        % Valid baseline inputs, reused across tests
        Nt
        inmk0
        radtheta
        epsf
        nmdz
        loc
    end

    methods (TestClassSetup)
        function addSrcToPath(~)
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));
        end
    end

    methods (TestMethodSetup)
        function buildBaselineInputs(testCase)
            Nt = 2;
            testCase.Nt       = Nt;
            testCase.inmk0    = 2*pi/600;
            testCase.radtheta = 0;

            ef = zeros(1, 4*Nt + 1);
            ef(1, 2*Nt + 1) = 2.25;
            testCase.epsf = ef;
            testCase.nmdz = 200;

            loc = DefaultLoc;
            loc.Nt   = Nt;
            loc.nmLx = 300;
            loc.pol  = 1;
            loc.nsa  = 1.0;
            testCase.loc = loc;
        end
    end

    % ======================================================================
    % RCWA
    % ======================================================================
    methods (Test)

        function testRCWARejectsNegativeWavenumber(testCase)
            testCase.verifyError( ...
                @() RCWA(-1, testCase.radtheta, testCase.epsf, testCase.nmdz, testCase.loc), ...
                'MATLAB:validators:mustBePositive');
        end

        function testRCWARejectsZeroWavenumber(testCase)
            testCase.verifyError( ...
                @() RCWA(0, testCase.radtheta, testCase.epsf, testCase.nmdz, testCase.loc), ...
                'MATLAB:validators:mustBePositive');
        end

        function testRCWARejectsComplexWavenumber(testCase)
            testCase.verifyError( ...
                @() RCWA(1i, testCase.radtheta, testCase.epsf, testCase.nmdz, testCase.loc), ...
                'MATLAB:validators:mustBeReal');
        end

        function testRCWARejectsNonScalarWavenumber(testCase)
            testCase.verifyError( ...
                @() RCWA([1, 2], testCase.radtheta, testCase.epsf, testCase.nmdz, testCase.loc), ...
                'MATLAB:validation:IncompatibleSize');
        end

        function testRCWARejectsNegativeLayerThickness(testCase)
            testCase.verifyError( ...
                @() RCWA(testCase.inmk0, testCase.radtheta, testCase.epsf, -100, testCase.loc), ...
                'MATLAB:validators:mustBePositive');
        end

        function testRCWARejectsEpsfNmdzSizeMismatch(testCase)
            % epsf has 1 row but nmdz has 2 elements
            testCase.verifyError( ...
                @() RCWA(testCase.inmk0, testCase.radtheta, testCase.epsf, [100, 200], testCase.loc), ...
                'RCWA:sizeMismatch');
        end

        % ======================================================================
        % BuildEps
        % ======================================================================

        function testBuildEpsRejectsNegativeWavelength(testCase)
            Nz      = 5;
            mat_cat = zeros(1, Nz);
            nmx     = linspace(-100, 100, 10);
            nmz     = linspace(0, 100, Nz);
            Eg      = zeros(Nz, 1);
            testCase.verifyError( ...
                @() BuildEps(mat_cat, nmx, nmz, -600, Eg, testCase.loc), ...
                'MATLAB:validators:mustBePositive');
        end

        function testBuildEpsRejectsNegativeMatCat(testCase)
            Nz      = 5;
            mat_cat = -1 * ones(1, Nz);  % invalid: negative category
            nmx     = linspace(-100, 100, 10);
            nmz     = linspace(0, 100, Nz);
            Eg      = zeros(Nz, 1);
            testCase.verifyError( ...
                @() BuildEps(mat_cat, nmx, nmz, 600, Eg, testCase.loc), ...
                'MATLAB:validators:mustBeNonnegative');
        end

        function testBuildEpsRejectsEgNmzSizeMismatch(testCase)
            Nz      = 5;
            mat_cat = zeros(1, Nz);
            nmx     = linspace(-100, 100, 10);
            nmz     = linspace(0, 100, Nz);
            Eg_bad  = zeros(Nz + 1, 1);  % one extra element
            testCase.verifyError( ...
                @() BuildEps(mat_cat, nmx, nmz, 600, Eg_bad, testCase.loc), ...
                'BuildEps:sizeMismatch');
        end

        % ======================================================================
        % BuildEpsF
        % ======================================================================

        function testBuildEpsFRejectsEpsxzNmzSizeMismatch(testCase)
            Nz    = 5;
            Nx    = 10;
            nmz   = linspace(0, 100, Nz);
            epsxz = ones(Nx, Nz + 1);  % wrong number of columns
            testCase.verifyError( ...
                @() BuildEpsF(nmz, epsxz, testCase.loc), ...
                'BuildEpsF:sizeMismatch');
        end

        function testBuildEpsFRejectsNonStructLoc(testCase)
            nmz   = linspace(0, 100, 5);
            epsxz = ones(10, 5);
            testCase.verifyError( ...
                @() BuildEpsF(nmz, epsxz, 42), ...  % scalar, not struct
                'MATLAB:validation:IncompatibleClass');
        end

        % ======================================================================
        % specular
        % ======================================================================

        function testSpecularRejectsNegativeWavenumber(testCase)
            Nmode = 2*testCase.Nt + 1;
            R = zeros(Nmode, 1);
            T = zeros(Nmode, 1);
            testCase.verifyError( ...
                @() specular(R, T, -1, testCase.radtheta, testCase.loc), ...
                'MATLAB:validators:mustBePositive');
        end

        function testSpecularRejectsRTLengthMismatch(testCase)
            Nmode = 2*testCase.Nt + 1;
            R = zeros(Nmode, 1);
            T = zeros(Nmode + 1, 1);  % one extra element
            testCase.verifyError( ...
                @() specular(R, T, testCase.inmk0, testCase.radtheta, testCase.loc), ...
                'specular:sizeMismatch');
        end

    end
end
