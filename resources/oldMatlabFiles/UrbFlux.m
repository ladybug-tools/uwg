function [UCM,UBL,BEM] = UrbFlux(UCM,UBL,BEM,forc,parameter,simTime,RSM)

    % Calculate the surface heat fluxes
    T_can = UCM.canTemp;
    Cp = parameter.cp;
    UCM.Q_roof = 0;
    sigma = 5.67e-8;        % Stephan-Boltzman constant
    UCM.roofTemp = 0;       % Average urban roof temperature
    UCM.wallTemp = 0;       % Average urban wall temperature
    
    for j = 1:numel(BEM)

        % Building energy model
        BEM(j).building = BEMCalc(BEM(j).building,UCM,BEM(j),forc,parameter,simTime);
        BEM(j).ElecTotal = BEM(j).building.ElecTotal * BEM(j).fl_area;

        % Update roof infra calc
        e_roof = BEM(j).roof.emissivity;
        T_roof = BEM(j).roof.layerTemp(1);     
        BEM(j).roof.infra = e_roof*(forc.infra - sigma*T_roof^4);
        
        % update wall infra calc (road done later)
        e_wall = BEM(j).wall.emissivity;
        T_wall = BEM(j).wall.layerTemp(1);
        [~,BEM(j).wall.infra]= InfraCalcs(UCM,forc, UCM.road.emissivity,e_wall,UCM.roadTemp,T_wall);
        
        % Update element temperatures
        BEM(j).mass.layerTemp = Conduction(BEM(j).mass,simTime.dt,BEM(j).building.fluxMass,1,0,BEM(j).building.fluxMass);
        BEM(j).roof = SurfFlux(BEM(j).roof,forc,parameter,simTime,UCM.canHum,T_can,max(forc.wind,UCM.canWind),1,BEM(j).building.fluxRoof);
        BEM(j).wall = SurfFlux(BEM(j).wall,forc,parameter,simTime,UCM.canHum,T_can,UCM.canWind,1,BEM(j).building.fluxWall);

        % Note the average wall & roof temperature
        UCM.wallTemp = UCM.wallTemp + BEM(j).frac*BEM(j).wall.layerTemp(1);
        UCM.roofTemp = UCM.roofTemp + BEM(j).frac*BEM(j).roof.layerTemp(1);

    end
    
    % Update road infra calc (assume walls have similar emissivity, so use the last one)
    [UCM.road.infra,~] = InfraCalcs(UCM,forc,UCM.road.emissivity,e_wall,UCM.roadTemp,UCM.wallTemp);
    UCM.road = SurfFlux(UCM.road,forc,parameter,simTime,UCM.canHum,T_can,UCM.canWind,2,0.);
    UCM.roadTemp = UCM.road.layerTemp(1);

    % Sensible & latent heat flux (total)
    UCM.latHeat = UCM.latHeat + UCM.latAnthrop + UCM.treeLatHeat + UCM.road.lat*(1-UCM.bldDensity);
 
    % ---------------------------------------------------------------------
    % Advective heat flux to UBL from VDM
    % ---------------------------------------------------------------------
    forDens = 0;
    intAdv1 = 0;
    intAdv2 = 0;
    for iz=1:RSM.nzfor
        forDens = forDens + RSM.densityProfC(iz)*RSM.dz(iz)/...
            (RSM.z(RSM.nzfor)+RSM.dz(RSM.nzfor)/2);
        intAdv1 = intAdv1 + RSM.windProf(iz)*RSM.tempProf(iz)*RSM.dz(iz);
        intAdv2 = intAdv2 + RSM.windProf(iz)*RSM.dz(iz);
    end
    UBL.advHeat = UBL.paralLength*Cp*forDens*(intAdv1-UBL.ublTemp*intAdv2)/UBL.urbArea;
    
    % ---------------------------------------------------------------------
    % Convective heat flux to UBL from UCM (see Appendix - Bueno (2014))
    % ---------------------------------------------------------------------
    zrUrb = 2*UCM.bldHeight;
    zref = RSM.z(RSM.nzref);    % Reference height

    % Reference wind speed & canyon air density
    windUrb = forc.wind*log(zref/RSM.z0r)/log(parameter.windHeight/RSM.z0r)*...
        log(zrUrb/UCM.z0u)/log(zref/UCM.z0u);
    dens = forc.pres/(1000*0.287042*T_can*(1.+1.607858*UCM.canHum));
    
    % Friction velocity
    UCM.ustar = parameter.vk*windUrb/log((zrUrb-UCM.l_disp)/UCM.z0u);
    
    % Convective scaling velocity
    wstar = (parameter.g*max(UCM.sensHeat,0)*zref/dens/Cp/T_can)^(1/3);
    UCM.ustarMod = max(UCM.ustar,wstar);        % Modified friction velocity
    UCM.uExch = parameter.exCoeff*UCM.ustarMod; % Exchange velocity

    % Canyon wind speed, Eq. 27 Chp. 3 Hanna and Britter, 2002
    % assuming CD = 1 and lambda_f = verToHor/4
    UCM.canWind = UCM.ustarMod*(UCM.verToHor/8)^(-1/2);

    % Canyon turbulent velocities
    UCM.turbU = 2.4*UCM.ustarMod;
    UCM.turbV = 1.9*UCM.ustarMod;
    UCM.turbW = 1.3*UCM.ustarMod;
    
    % Urban wind profile
    for iz=1:RSM.nzref
        UCM.windProf(iz) = UCM.ustar/parameter.vk*...
            log((RSM.z(iz)+UCM.bldHeight-UCM.l_disp)/UCM.z0u);
    end
end