classdef Param
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % system parameters
        dayBLHeight;    % daytime mixing height, orig = 700
        nightBLHeight;  % Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
        refHeight;      % Reference height at which the vertical profile 
                        %   of potential temperature is vertical
        tempHeight;     % Temperature measuremnt height at the weather station (m)
        windHeight;     % Air velocity measuremnt height at the weather station (m)
        circCoeff;      % Wind scaling coefficient
        dayThreshold;   % heat flux threshold for daytime conditions (W m-2)
        nightThreshold; % heat flux threshold for nighttime conditions (W m-2)
        treeFLat;       % latent fraction of trees
        grassFLat;      % latent fraction of grass
        vegAlbedo;      % albedo of vegetation
        vegStart;       % begin month for vegetation participation
        vegEnd;         % end month for vegetation participation
        nightSetStart;  % begin hour for night thermal set point schedule
        nightSetEnd;    % end hour for night thermal set point schedule
        windMin;        % minimum wind speed (m s-1)
        wgmax;          % maximum film water depth on horizontal surfaces (m)
        exCoeff;        % exchange velocity coefficient
        maxdx;          % maximum discretization length for the UBL model (m)
        % physical parameters
        g;              % gravity
        cp;             % heat capacity for air ((J/kg.K))
        vk;             % von karman constant
        r;              % gas constant
        rv;
        lv;             % latent heat of evaporation
        pi;             % pi
        sigma ;         % Stefan Boltzmann constant
        waterDens;      % water density
        lvtt;
        tt;
        estt;
        cl;
        cpv;
        b;              % Coefficients derived by Louis (1979)
        cm; 
        colburn;        % (Pr/Sc)^(2/3) for Colburn analogy in water evaporation
    end

    methods
        function obj = Param(dayBLHeight, nightBLHeight, refHeight, tempHeight, windHeight, ...
                circCoeff, dayThreshold, nightThreshold, treeFLat, grassFLat, vegAlbedo, vegStart, ...
                vegEnd, nightSetStart, nightSetEnd, windMin, wgmax, exCoeff, maxdx, g, ...
                cp, vk, r, rv, lv, pi, sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn)
            % class constructor
            if(nargin > 0)
                obj.dayBLHeight = dayBLHeight;
                obj.nightBLHeight = nightBLHeight;
                obj.refHeight = refHeight;
                obj.tempHeight = tempHeight; 
                obj.windHeight = windHeight;
                obj.circCoeff =  circCoeff;
                obj.dayThreshold = dayThreshold;
                obj.nightThreshold = nightThreshold; 
                obj.treeFLat = treeFLat; 
                obj.grassFLat = grassFLat;
                obj.vegAlbedo = vegAlbedo;
                obj.vegStart = vegStart;
                obj.vegEnd = vegEnd;
                obj.nightSetStart = nightSetStart;
                obj.nightSetEnd = nightSetEnd;
                obj.windMin = windMin;
                obj.wgmax = wgmax;
                obj.exCoeff = exCoeff;
                obj.maxdx = maxdx;
                obj.g = g;
                obj.cp = cp;
                obj.vk = vk;
                obj.r = r;
                obj.rv = rv;
                obj.lv = lv;
                obj.pi = pi;
                obj.sigma = sigma;
                obj.waterDens = waterDens;
                obj.lvtt = lvtt;
                obj.tt = tt;
                obj.estt = estt;
                obj.cl = cl;
                obj.cpv = cpv;
                obj.b   = b;
                obj.cm  = cm;
                obj.colburn = colburn;
            end
        end          
    end
end

