classdef GratingTest < matlab.unittest.TestCase
    methods (Test)
        function squareReliefMatchesLegacy(tc)
            cfg = rcwa.Config('type', 0, 'nmLx', 400, 'zeta', 0.5, 'nmLg', 50);
            sq  = rcwa.gratings.Square('nmLx', 400, 'zeta', 0.5, 'nmLg', 50);
            nmx = linspace(-200, 200, 51);
            tc.verifyEqual(sq.relief(nmx), gx(nmx, cfg.toStruct()), ...
                'AbsTol', 1e-12);
        end

        function sinusoidalReliefMatchesLegacy(tc)
            cfg = rcwa.Config('type', 1, 'nmLx', 400, 'zeta', 0.6, 'nmLg', 75);
            si  = rcwa.gratings.Sinusoidal('nmLx', 400, 'zeta', 0.6, 'nmLg', 75);
            nmx = linspace(-200, 200, 101);
            tc.verifyEqual(si.relief(nmx), gx(nmx, cfg.toStruct()), ...
                'AbsTol', 1e-12);
        end

        function pyramidReliefMatchesLegacy(tc)
            cfg = rcwa.Config('type', 2, 'nmLx', 400, 'zeta', 0.5, 'nmLg', 30);
            py  = rcwa.gratings.Pyramid('nmLx', 400, 'zeta', 0.5, 'nmLg', 30);
            nmx = linspace(-200, 200, 51);
            tc.verifyEqual(py.relief(nmx), gx(nmx, cfg.toStruct()), ...
                'AbsTol', 1e-12);
        end

        function squarePermittivityCommutes(tc)
            % Square + uniform z within grating region: epsAt returns
            % epsm where relief covers, epsd elsewhere.
            sq = rcwa.gratings.Square('nmLx', 400, 'zeta', 0.5, 'nmLg', 50);
            nmx = linspace(-200, 200, 101);
            slice = sq.epsAt(nmx, 110, 100, 1+2i, 5+0.1i, 0);
            % at z=110 (inside the relief which is 50 above mirror at 100)
            % under the grating ridge (|x|<100) we should be in metal:
            tc.verifyEqual(slice(abs(nmx) < 99), ...
                (1+2i) * ones(1, sum(abs(nmx) < 99)), 'AbsTol', 1e-12);
            % outside, dielectric:
            tc.verifyEqual(slice(abs(nmx) > 101), ...
                (5+0.1i) * ones(1, sum(abs(nmx) > 101)), 'AbsTol', 1e-12);
        end

        function factoryReturnsRightClass(tc)
            names = {'rcwa.gratings.Square', ...
                     'rcwa.gratings.Sinusoidal', ...
                     'rcwa.gratings.Pyramid', ...
                     'rcwa.gratings.Spherical'};
            for typ = 0:3
                cfg = rcwa.Config('type', typ);
                g   = rcwa.gratings.Grating.fromConfig(cfg);
                tc.verifyClass(g, names{typ+1});
            end
        end

        function squareFourierExplicitClosedForm(tc)
            % The exact Fourier coefficients of a rectangular profile of
            % duty cycle zeta and contrast (epsm-epsd) are known:
            %    f_n = (epsm-epsd) * sin(n*pi*zeta) / (n*pi)   for n != 0
            sq = rcwa.gratings.Square('zeta', 0.5);
            n = [-3 -2 -1 1 2 3];
            f = sq.fourierExplicit(n, 1+0i, 5+0i);
            expected = (1 - 5) .* sin(n * pi * 0.5) ./ (n * pi);
            tc.verifyEqual(f, expected, 'AbsTol', 1e-12);
        end
    end
end
