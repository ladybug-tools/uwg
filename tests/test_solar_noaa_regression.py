"""NOAA equation checks at the supplied instant, without changing its timestamp."""
import datetime
import math
from types import SimpleNamespace

import pytest
from uwg import SimParam, SolarCalcs


@pytest.mark.parametrize('site', [(1.35, 103.8, 8), (49.4, 8.7, 1),
                                 (39.74, -105, -7), (-33.86, 151.2, 10)])
@pytest.mark.parametrize('month,day', [(3, 21), (6, 21), (9, 22), (12, 21)])
@pytest.mark.parametrize('hour', [0, 6, 9, 12, 15, 18, 23])
def test_noaa_geometric_zenith_at_local_standard_time(site, month, day, hour):
    # NOAA General Solar Position Calculations, gml.noaa.gov/grad/solcalc/solareqns.PDF
    lat, lon, timezone = site
    time = SimParam(300, 3600, month, day, 1)
    time.secDay = hour * 3600
    solar = SolarCalcs(SimpleNamespace(canAspect=1), [], time,
                       SimpleNamespace(lat=lat, lon=lon, gmt=timezone), None, None, None)
    solar.solarangles()
    doy = datetime.date(2025, month, day).timetuple().tm_yday
    gamma = 2 * math.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eq = 229.18 * (.000075 + .001868 * math.cos(gamma) - .032077 * math.sin(gamma)
                   - .014615 * math.cos(2 * gamma) - .040849 * math.sin(2 * gamma))
    dec = (.006918 - .399912 * math.cos(gamma) + .070257 * math.sin(gamma)
           - .006758 * math.cos(2 * gamma) + .000907 * math.sin(2 * gamma)
           - .002697 * math.cos(3 * gamma) + .00148 * math.sin(3 * gamma))
    hour_angle = math.radians((hour * 60 + eq + 4 * lon - 60 * timezone) / 4 - 180)
    expected = math.acos(math.sin(math.radians(lat)) * math.sin(dec) +
                         math.cos(math.radians(lat)) * math.cos(dec) * math.cos(hour_angle))
    assert solar.ad == pytest.approx(gamma, abs=1e-12)
    assert math.degrees(solar.zenith) == pytest.approx(math.degrees(expected), abs=.001)
    assert time.secDay == hour * 3600
