"""Functions for calculating canyon infrared radiation."""

SIGMA = 5.67e-8  # Stephen-Boltzman const


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
    # road_wall_conf: configuration factor (view factor) for road to wall
    road_wall_conf = (1. - UCM.roadConf)
    # wall_road_conf: wall to road VF same as wall-sky configuration factors
    # (sky view factor)
    wall_road_conf = UCM.wallConf

    # Calculate radiation of unshaded road, accounting for radiation exchange from wall
    _road_rad = \
        e_road * UCM.roadConf * (1. - UCM.roadShad) * (forc.infra - SIGMA * T_road ** 4.)
    _wall_to_road_rad = \
        (1. - UCM.roadShad) * e_wall * e_road * SIGMA * road_wall_conf * (T_wall ** 4. - T_road ** 4.)

    infra_road = _road_rad + _wall_to_road_rad

    # Calculate radiation of wall, accounting for radiation exchange from unshaded road
    _road_rad = e_wall * UCM.wallConf * (forc.infra - SIGMA * T_wall ** 4.)
    _road_to_wall_rad = \
        (1. - UCM.roadShad) * e_wall * e_road * SIGMA * wall_road_conf * (T_road ** 4. - T_wall ** 4.)

    infra_wall = _road_rad + _road_to_wall_rad

    return infra_road, infra_wall
