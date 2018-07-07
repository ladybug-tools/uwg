def infracalcs(UCM,forc,e_road,e_wall,T_road,T_wall):
    # Calculate IR surface flux
    sigma = 5.67e-8                         # Stephen-Boltzman const
    road_wall_conf = (1. - UCM.roadConf)    # configuration factor (view factor) for road to wall
    wall_road_conf = UCM.wallConf           # wall to road VF same as wall-sky configuration factors (sky view factor)

    # Calculate radiation of unshaded road, accounting for radiation exchange from wall
    infra_road = e_road * UCM.roadConf * (1.-UCM.roadShad) * (forc.infra - sigma*T_road**4.) + (1.-UCM.roadShad) * e_wall * e_road * sigma * road_wall_conf * (T_wall**4.-T_road**4.)
    # Calculate radiation of wall, accounting for radiation exchange from unshaded road
    infra_wall = e_wall * UCM.wallConf * (forc.infra - sigma*T_wall**4.) + (1.-UCM.roadShad) * e_wall * e_road * sigma * wall_road_conf * (T_road**4.-T_wall**4.)

    return infra_road, infra_wall
