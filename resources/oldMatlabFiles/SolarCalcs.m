function [rural,UCM,BEM] = SolarCalcs(UCM,BEM,simTime,RSM,forc,parameter,rural)

    dir = forc.dir;     % Direct sunlight (perpendicular to the sun's ray)
    dif = forc.dif;     % Diffuse sunlight
    
    if dir + dif > 0
        
        % Solar angles
        [zenith, tanzen, critOrient] = SolarAngles(UCM.canAspect,simTime,RSM.lon,RSM.lat,RSM.GMT);
        horSol = max(cos(zenith)*dir,0);            % Direct horizontal radiation

        % Fractional terms for wall & road
        Kw_term = min(abs(1/UCM.canAspect*(0.5-critOrient/pi)+1/pi*tanzen*(1-cos(critOrient))),1);
        Kr_term = min(abs(2*critOrient/pi-(2/pi*UCM.canAspect*tanzen)*(1-cos(critOrient))),1-2*UCM.canAspect*Kw_term);

        % Direct and diffuse solar radiation
        bldSol = horSol*Kw_term+UCM.wallConf*dif;   % Assume trees are shorter than buildings
        roadSol = horSol*Kr_term+UCM.roadConf*dif;

        % Solar reflections
        if simTime.month < parameter.vegStart || simTime.month > parameter.vegEnd
            alb_road = UCM.road.albedo;
        else
            alb_road = UCM.road.albedo*(1-UCM.road.vegCoverage)+parameter.vegAlbedo*UCM.road.vegCoverage;
        end
        
        % First set of reflections
        rr = alb_road*roadSol;
        rw = UCM.alb_wall*bldSol;
        
        % bounces
        mr = (rr+(1-UCM.roadConf)*alb_road*...
            (rw+UCM.wallConf*UCM.alb_wall*rr))/...
            (1-(1-2*UCM.wallConf)*UCM.alb_wall+...
            (1-UCM.roadConf)*UCM.wallConf*alb_road*UCM.alb_wall);
        mw = (rw+UCM.wallConf*UCM.alb_wall*rr)/...
            (1-(1-2*UCM.wallConf)*UCM.alb_wall+...
            (1-UCM.roadConf)*UCM.wallConf*alb_road*UCM.alb_wall);

        % Receiving solar, including bounces
        UCM.road.solRec = roadSol+(1-UCM.roadConf)*mw;
        for j = 1:numel(BEM)
            BEM(j).roof.solRec = horSol + dif;  
            BEM(j).wall.solRec = bldSol+(1-2*UCM.wallConf)*mw+UCM.wallConf*mr; 
        end
        
        rural.solRec = horSol + dif;            % Solar received by rural 
        UCM.SolRecRoof = horSol + dif;          % Solar received by roof
        UCM.SolRecRoad = UCM.road.solRec;       % Solar received by road
        UCM.SolRecWall = bldSol+(1-2*UCM.wallConf)*UCM.road.albedo*roadSol;    % Solar received by wall    

        % Vegetation heat (per m^2 of veg)
        UCM.treeSensHeat = (1-parameter.vegAlbedo)*(1-parameter.treeFLat)*UCM.SolRecRoad;
        UCM.treeLatHeat = (1-parameter.vegAlbedo)*parameter.treeFLat*UCM.SolRecRoad;
        
    else    % No Sun
        
        UCM.road.solRec = 0;
        rural.solRec = 0;

        for j = 1:numel(BEM)
            BEM(j).roof.solRec = 0;
            BEM(j).wall.solRec = 0;
        end
        UCM.SolRecRoad = 0;         % Solar received by road
        UCM.SolRecRoof = 0;         % Solar received by roof
        UCM.SolRecWall = 0;         % Solar received by wall    
        UCM.treeSensHeat = 0;
        UCM.treeLatHeat = 0;
        
    end
end

function [zenith, tanzen, theta0] = SolarAngles (canAspect,simTime,lon,lat,GMT)

    % Calculation based on NOAA 
    month = simTime.month;
    day = simTime.day;
    secDay = simTime.secDay;
    inobis = simTime.inobis;
    ut = mod(24.0+mod(secDay/3600.,24.0),24.0);
    
    ibis = inobis;
    for JI=2:12
      ibis(JI) = inobis(JI)+1;
    end
    
    date = day + inobis(month) - 1;         % Julian day of the year
    ad = 2.0*pi/365*(date-1+(ut-12/24));    % Fractional year (rad)
    
    eqtime = 229.18*(0.000075+0.001868*cos(ad)-0.032077*sin(ad) - ...
        0.01461*cos(2*ad)-0.040849*sin(2*ad));

    % Declination angle
    decsol = 0.006918-0.399912*cos(ad)+0.070257*sin(ad)...
             -0.006758*cos(2.*ad)+0.000907*sin(2.*ad)...
             -0.002697*cos(3.*ad)+0.00148 *sin(3.*ad);
    
    time_offset = eqtime - 4*lon + 60*(GMT);
    tst = secDay + time_offset*60;
    ha = (tst/4/60-180)*pi/180;
    zlat = lat*(pi/180.);   % change angle units

    % Calculate zenith solar angle
    zenith = acos(sin(zlat)*sin(decsol) + cos(zlat)*cos(decsol)*cos(ha));
    
    % tangente of solar zenithal angle
    if (abs(0.5*pi-zenith) <  1.E-6)
        if(0.5*pi-zenith > 0.)  
            tanzen = tan(0.5*pi-1.E-6); 
        end
        if(0.5*pi-zenith <= 0.) 
            tanzen = tan(0.5*pi+1.E-6); 
        end
    elseif (abs(zenith) <  1.E-6)
        tanzen = sign(1.,zenith)*tan(1.E-6);
    else
        tanzen = tan(zenith);
    end
    
    % critical canyon angle for which solar radiation reaches the road
    theta0 = asin(min(abs( 1./tanzen)/canAspect, 1. ));
end