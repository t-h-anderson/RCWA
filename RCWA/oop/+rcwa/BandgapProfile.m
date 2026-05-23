classdef BandgapProfile
    % Bandgap Eg(z) across the p-i-n junction.
    %
    % Replaces EgProfile.m. The i-layer can carry a periodic perturbation:
    %     Eg(z) = Eg0 + A * (0.5*(sin(2*pi*((z-z_i)*kappa/Li + phi)) + 1))^alpha
    % The p- and n-layers default to a fixed wide-bandgap "window" value of
    % 1.95 eV (faithful to the original).
    %
    % NOTE on BUG-001: the original EgProfile uses (z - nmLn) inside the
    % sinusoid, which centres the perturbation on the n-layer thickness
    % instead of the i-layer start (= nmLp). This class exposes a
    % `legacy` flag that reproduces the bug for thesis-reproducibility.

    properties
        Eg0    (1,1) double = 1.6
        A      (1,1) double = 0
        kappa  (1,1) double = 0
        phi    (1,1) double = 0
        alpha  (1,1) double = 0
        EgPN   (1,1) double = 1.95   % bandgap in p- and n-layers
        nmLp   (1,1) double = 15
        nmLi   (1,1) double = 200
        nmLn   (1,1) double = 15
        legacy (1,1) logical = true  % reproduce BUG-001 by default
    end

    methods
        function obj = BandgapProfile(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function obj = fromConfig_(obj, cfg)
            % Internal: pull profile params from a Config.
            obj.Eg0  = cfg.Eg0;   obj.A     = cfg.A;
            obj.kappa = cfg.kappa; obj.phi  = cfg.phi;
            obj.alpha = cfg.alpha;
            obj.nmLp = cfg.nmLp;  obj.nmLi = cfg.nmLi; obj.nmLn = cfg.nmLn;
        end

        function Eg = evaluate(obj, nmz)
            % Eg(z) sampled at the row vector nmz.
            Eg = zeros(1, numel(nmz));

            player = sign(nmz - obj.nmLp);
            player = 0.5*(-player + abs(player));

            ilayer = sign(nmz - (obj.nmLp + obj.nmLi)) + player;
            ilayer = 0.5*(-ilayer + abs(ilayer));

            nlayer = sign(nmz - (obj.nmLp + obj.nmLi + obj.nmLn)) + ilayer + player;
            nlayer = 0.5*(-nlayer + abs(nlayer));

            Eg(nlayer == 1) = obj.EgPN;
            Eg(player == 1) = obj.EgPN;

            if obj.legacy
                offset = obj.nmLn;   % BUG-001 reproduction
            else
                offset = obj.nmLp;
            end
            zi = nmz(ilayer == 1);
            inner = 0.5 * (sin( 2*pi * (zi - offset) * obj.kappa / obj.nmLi ...
                                + 2*pi*obj.phi ) + 1);
            Eg(ilayer == 1) = obj.Eg0 + obj.A .* inner.^obj.alpha;
        end
    end

    methods (Static)
        function obj = fromConfig(cfg)
            obj = rcwa.BandgapProfile();
            obj = obj.fromConfig_(cfg);
        end
    end
end
