function [infra_road,infra_wall] = InfraCalcs(UCM,forc,e_road,e_wall,T_road,T_wall)
    % Calculate IR surface flux 
    sigma = 5.67e-8;                % Stephen-Boltzman const
    roadWallConf =(1 - UCM.roadConf);   
    wallRoadConf = UCM.wallConf;    % wall to road VF same as wall to sky
    
    infra_road = e_road*UCM.roadConf*(1-UCM.roadShad)*(forc.infra - sigma*T_road^4)+(1-UCM.roadShad)*e_wall*e_road*sigma*roadWallConf*(T_wall^4-T_road^4);
    infra_wall = e_wall*UCM.wallConf*(forc.infra - sigma*T_wall^4)+(1-UCM.roadShad)*e_wall*e_road*sigma*wallRoadConf*(T_road^4-T_wall^4);

end