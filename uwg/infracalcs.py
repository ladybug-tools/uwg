"""Functions for calculating canyon infrared radiation."""


def infracalcs(UCM, forc, e_road, e_wall, T_road, T_wall):
    """Calculate infrared radiation surface flux on canyon wall and road.

    Args:
        UCM: UCMDef object.
        forc: Forcing object.
        e_road: Number for road emissivity.
        e_wall: Number for wall emissivity.
        T_road: Number for road surface temperature (K).
        T_wall: Number for wall surface temperature (K).

    Returns:
        Tuple of two values

        - infra_road: Road infrared radiation (W)

        - infra_wall: Wall infrared radiation (W)
    """
    # sigma: Stephen-Boltzman const
    sigma = 5.67e-8
    # road_wall_conf: configuration factor (view factor) for road to wall
    road_wall_conf = (1. - UCM.roadConf)
    # wall_road_conf: wall to road VF same as wall-sky configuration factors
    # (sky view factor)
    wall_road_conf = UCM.wallConf

    # Calculate radiation of unshaded road, accounting for radiation exchange from wall
    infra_road = e_road * UCM.roadConf * (1. - UCM.roadShad) * \
        (forc.infra - sigma * T_road ** 4.) + (1. - UCM.roadShad) * \
        e_wall * e_road * sigma * road_wall_conf * (T_wall ** 4. - T_road ** 4.)

    # Calculate radiation of wall, accounting for radiation exchange from unshaded road
    infra_wall = e_wall * UCM.wallConf * (forc.infra - sigma * T_wall ** 4.) + \
        (1. - UCM.roadShad) * e_wall * e_road * sigma * wall_road_conf * \
        (T_road ** 4.0 - T_wall ** 4.0)

    return infra_road, infra_wall
