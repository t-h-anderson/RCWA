classdef InverseFourier
    % Reconstruction of a real-space field from its Fourier coefficients.
    %
    % Replaces InverseFourier.m. Two methods are available:
    %   method = 0 -> explicit Fourier sum sum_j f_j exp(i*k_j*x)
    %   method = 1 -> ifft + interpft path

    methods (Static)
        function y = reconstruct(x, fcomp, Lx, method, NTcap)
            fcomp_temp = squeeze(fcomp);
            terms = numel(fcomp_temp);
            NT = (terms - 1) / 2;

            if nargin >= 5 && ~isempty(NTcap)
                NT0 = min(NT, floor(NTcap));
            else
                NT0 = NT;
            end

            if method == 1
                fcomp_temp(2:2:end) = -fcomp_temp(2:2:end);
                FFTorder = (2*NT+1) * [fcomp_temp(NT+1:NT+1+NT0); ...
                                       fcomp_temp(NT-NT0+1:NT)];
                invff = ifft(FFTorder);
                y = interpft(invff, numel(x));
            else
                y = zeros(size(x));
                for j = -NT0:NT0
                    k = 2*pi * j / Lx;
                    y = y + fcomp(j + NT + 1) .* exp(1i * k * x);
                end
            end
        end
    end
end
