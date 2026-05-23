classdef InverseFourierTest < matlab.unittest.TestCase
    methods (Test)
        function explicitMatchesAnalytic(tc)
            % If only the n=0 mode is non-zero, the reconstruction should
            % be a uniform field equal to that mode.
            Lx = 400;
            Nt = 3;
            f  = zeros(2*Nt + 1, 1);
            f(Nt+1) = 2.5;     % zeroth mode
            x  = linspace(-Lx/2, Lx/2, 11);
            y  = rcwa.InverseFourier.reconstruct(x, f, Lx, 0);
            tc.verifyEqual(y, 2.5 * ones(size(x)), 'AbsTol', 1e-12);
        end

        function singlePositiveMode(tc)
            Lx = 400;
            Nt = 2;
            f  = zeros(2*Nt+1, 1);
            f(Nt+2) = 1;       % f_{+1} = 1
            x = linspace(-Lx/2, Lx/2, 100);
            y = rcwa.InverseFourier.reconstruct(x, f, Lx, 0);
            expected = exp(1i * 2*pi/Lx * x);
            tc.verifyEqual(y, expected, 'AbsTol', 1e-12);
        end

        function matchesLegacyInverseFourier(tc)
            % Cross-check the OOP reconstruction against the legacy
            % InverseFourier.m result for a random spectrum.
            rng(42);
            Lx = 250;
            Nt = 4;
            f  = randn(2*Nt+1, 1) + 1i * randn(2*Nt+1, 1);
            x  = linspace(-Lx/2, Lx/2, 33);
            y_oop = rcwa.InverseFourier.reconstruct(x, f, Lx, 0, Nt);
            y_leg = InverseFourier(x, f, Lx, 0, Nt);
            tc.verifyEqual(y_oop(:), y_leg(:), 'AbsTol', 1e-10);
        end
    end
end
