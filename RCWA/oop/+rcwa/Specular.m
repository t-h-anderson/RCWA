classdef Specular
    % Reflection / transmission specular and diffuse components.
    %
    % Replaces specular.m. Pure transforms, no I/O.

    methods (Static)
        function [rp0, rp1, tp0, tp1] = compute(cfg, R, T, inmk0, radtheta)
            Nt = cfg.Nt;
            kxn = cfg.nsa .* inmk0 .* sin(radtheta) ...
                  + (-Nt:Nt) * 2*pi / cfg.nmLx;
            kzn = sqrt(cfg.nsa^2 .* inmk0^2 - kxn.^2);

            denom = cfg.nsa * inmk0 * cos(radtheta);
            kzr   = real(kzn) / denom;

            RP2 = abs(R).^2;
            TP2 = abs(T).^2;

            RPp = RP2(:) .* kzr(:);   % element-wise, was diag-mult in original
            TPp = TP2(:) .* kzr(:);

            rp0 = RPp(Nt+1);
            rp1 = sum(RPp) - RPp(Nt+1);
            tp0 = TPp(Nt+1);
            tp1 = sum(TPp) - TPp(Nt+1);
        end
    end
end
