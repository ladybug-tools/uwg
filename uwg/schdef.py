"""Class for Schedule Definition."""


class SchDef(object):
    """Schedule class

    Args:
        elec: Matrix of numbers for weekly schedule for electricity (WD, Sat, Sun).
            (Default: None).
        gas: Matrix of numbers for weekly schedule for gas (WD, Sat, Sun).
            (Default: None)
        light: Matrix of numbers for weekly schedule for light (WD, Sat, Sun).
            (Default: None).
        occ: Matrix of numbers for weekly schedule for occupant (WD, Sat, Sun).
            (Default: None).
        cool: Matrix of numbers for weekly temperature schedule for cooling
            (WD, Sat, Sun). (Default: None).
        heat: Matrix of numbers for weekly temperature schedule for heating
            (WD, Sat, Sun). (Default: None).
        swh: Matrix of numbers for weekly hot water schedule (WD, Sat, Sun).
            (Default: None).
        bldtype: Number between 0 and 15 corresponding to the following building
            types: FullServiceRestaurant (0), Hospital (1), LargeHotel (2),
            LargeOffice (3), MediumOffice (4), MidRiseApartment (5), OutPatient (6),
            PrimarySchool (7), QuickServiceRestaurant (8), SecondarySchool (9),
            SmallHotel (10), SmallOffice (11), StandaloneRetail (12), StripMall (13),
            SuperMarket (14), Warehouse (15). Additional building types can be defined
            with a number greater then 15. Default is None.
        builtera: Number between 0 and 2 corresponding to the following built eras:
            Pre-1980s (0), Post1980s (1), New construction (2). Default is None.
        zonetype: Number between 0 and 15 corresponding to the following zone types:
            1A-Miami (0), 2A-Houston (1), 2B-Phoenix (2), 3A-Atlanta (3),
            3B-CA-Los Angeles (4), 3B-Las Vegas (5), 3C (San Francisco) (6),
            4A (Baltimore) (7), 4B (Albuquerque) (8), 4C (Seattle) (9),
            5A (Chicago) (10), 5B (Boulder) (11), 6A (Minneapolis) (12),
            6B (Helena) (13), 7 (Duluth) (14), 8 (Fairbanks) (15). Default is None.

    Properties:
        * elec
        * gas
        * light
        * occ
        * cool
        * heat
        * swh
        * bldtype
        * builtera
        * zonetype
    """

    def __init__(self, elec, gas, light, occ, cool, heat, swh, bldtype=None,
                 builtera=None, zonetype=None):
        self.elec = SchDef.check_week_validity(elec, 'elec')
        self.gas = SchDef.check_week_validity(gas, 'gas')
        self.light = SchDef.check_week_validity(light, 'light')
        self.occ = SchDef.check_week_validity(occ, 'occ')
        self.cool = SchDef.check_week_validity(cool, 'cool')
        self.heat = SchDef.check_week_validity(heat, 'heat')
        self.swh = SchDef.check_week_validity(swh, 'swh')

        # Properties to be set in readDOE
        self.bldtype = bldtype  # DOE reference building type
        self.builtera = builtera  # pre80, pst80, new
        self.zonetype = zonetype  # climate zone number

    @classmethod
    def from_dict(cls, data):
        """Create a SchDef object from a dictionary.

        Args:
            data: A SchDef dictionary following the format below.

        .. code-block:: python

            _example_week = [
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6,
                 0.7, 0.8, 0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                 0.6, 0.7, 0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday

            {
            "elec": _example_week,  # 3 x 24 matrix of hourly values
            "gas": _example_week,  # 3 x 24 matrix of hourly values
            "light": _example_week,  # 3 x 24 matrix of hourly values
            "occ": _example_week,  # 3 x 24 matrix of hourly values
            "cool": _example_week,  # 3 x 24 matrix of hourly values
            "heat": _example_week,  # 3 x 24 matrix of hourly values
            "swh": _example_week,  # 3 x 24 matrix of hourly values
            "bldtype": 0,  # building type index
            "builtera": 1,  # built era index
            "zonetype": 15  # zone type index
            }
        """
        assert data['type'] == 'SchDef', 'Expected ' \
            'SchDef dictionary. Got {}.'.format(data['type'])

        sch = cls(elec=data['elec'], gas=data['gas'], light=data['light'],
                  occ=data['occ'], cool=data['cool'], heat=data['heat'],
                  swh=data['swh'])

        sch.bldtype = data['bldtype']
        sch.zonetype = data['zonetype']
        sch.builtera = data['builtera']

    def to_dict(self):
        """SchDef dictionary representation."""
        base = {'type': 'SchDef'}
        base['elec'] = self.elec
        base['gas'] = self.gas
        base['light'] = self.light
        base['occ'] = self.occ
        base['cool'] = self.cool
        base['heat'] = self.heat
        base['swh'] = self.swh
        base['bldtype'] = self.bldtype
        base['zonetype'] = self.zonetype
        base['builtera'] = self.builtera
        return base

    @staticmethod
    def check_week_validity(week, name):
        assert isinstance(week, (list, tuple)), 'The {} property must be a ' \
            'list or tuple. Got {}.'.format(name, week)

        assert len(week) == 3, 'The {} property must be a 3 x 24 matrix. Got ' \
            '{} rows.'.format(name, len(week))

        for i, day in enumerate(week):
            assert len(day) == 24, 'The {} property must be a 3 x 24 ' \
                'matrix. Got {} columns for row {}.'.format(name, len(day), i)
            for val in day:
                assert isinstance(val, (float, int)), 'The {} property ' \
                    'must be a 3 x 23 matrix of numbers. Got: {}.'.format(name, day)
        return week

    def __repr__(self):
        return 'Schedule, bldtype: {}\n zonetype: {}\n builtera: {}\n '\
            'Setpoint schedules:(weekday from 8:00 - 18:00)\n Heating: {}\n ' \
            'Cooling: {}'.format(
                self.bldtype, self.zonetype, self.builtera,
                'Null' if self.heat is None else self.heat[0][7:17],
                'Null' if self.cool is None else self.cool[0][7:17])
