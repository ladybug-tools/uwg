"""Class for Schedule Definition."""


class SchDef(object):
    """Schedule definition class.

    Args:
        elec: Weekly schedule of fractional electricity plug process load. Weekly
            schedule consists of three lists of 24 values representing hours in
            weekday, Saturday, and Sunday.
        gas: Matrix of numbers for weekly schedule for gas (WD, Sat, Sun).
        light: Matrix of numbers for weekly schedule for light (WD, Sat, Sun).
        occ: Matrix of numbers for weekly schedule for occupant (WD, Sat, Sun).
        cool: Matrix of numbers for weekly temperature schedule for cooling
            (WD, Sat, Sun).
        heat: Matrix of numbers for weekly temperature schedule for heating
            (WD, Sat, Sun).
        swh: Matrix of numbers for weekly hot water schedule (WD, Sat, Sun).
        q_elec: Number for maximum electrical plug process load [W/m2].
        q_gas: Number for maximum gas process load per unit area [W/m2].
        q_light: Number for maximum light process load per unit area [W/m2].
        n_occ: Number for maximum number of occupants per unit area [person/m2].
        vent: Number for maximum ventilation rate per unit area [m3/s/m2].
        v_swh: Number for maximum hot water rate per unit area [L/hr/m2].
        bldtype: Number between 0 and 15 corresponding to the following building
            types: FullServiceRestaurant (0), Hospital (1), LargeHotel (2),
            LargeOffice (3), MediumOffice (4), MidRiseApartment (5), OutPatient (6),
            PrimarySchool (7), QuickServiceRestaurant (8), SecondarySchool (9),
            SmallHotel (10), SmallOffice (11), StandaloneRetail (12), StripMall (13),
            SuperMarket (14), Warehouse (15). Additional building types can be defined
            with a number greater then 15. This value is used to reference the fraction
            of urban area the SchDef object defines in the UWG bld matrix.
        builtera: Number between 0 and 2 corresponding to the following built eras:
            Pre-1980s (0), Post1980s (1), New construction (2). This value is used to
            reference the fraction of urban area the SchDef object defines in the UWG
            bld matrix.

    Properties:
        * elec
        * gas
        * light
        * occ
        * cool
        * heat
        * swh
        * q_elec
        * q_gas
        * q_light
        * n_occ
        * vent
        * v_swh
        * bldtype
        * builtera
        * zonetype
    """

    def __init__(self, elec, gas, light, occ, cool, heat, swh, q_elec, q_gas, q_light,
                 n_occ, vent, v_swh, bldtype, builtera):
        self.elec = SchDef.check_week_validity(elec, 'elec')
        self.gas = SchDef.check_week_validity(gas, 'gas')
        self.light = SchDef.check_week_validity(light, 'light')
        self.occ = SchDef.check_week_validity(occ, 'occ')
        self.cool = SchDef.check_week_validity(cool, 'cool')
        self.heat = SchDef.check_week_validity(heat, 'heat')
        self.swh = SchDef.check_week_validity(swh, 'swh')
        self.q_elec = q_elec
        self.q_gas = q_gas
        self.q_light = q_light
        self.n_occ = n_occ
        self.vent = vent
        self.v_swh = v_swh
        # Properties to be set in readDOE
        self.bldtype = bldtype  # DOE reference building type
        self.builtera = builtera  # pre80, pst80, new
        self.zonetype = None  # climate zone number (only used in testing).

    # @property
    # def elec(self):
    #     """Get or set weekly schedule of fractional electricity plug process load.

    #     Weekly schedule consists of three lists of 24 values representing hours in
    #     weekday, Saturday, and Sunday.
    #     """
    #     return self._elec

    # @elec.setter
    # def elec(self, value):
    #     self._elec = SchDef.check_week_validity(value, 'elec')

    # @property
    # def gas(self):
    #     """Get or set weekly schedule of fractional gas process load.

    #     Weekly schedule consists of three lists of 24 values representing hours in
    #     weekday, Saturday, and Sunday.
    #     """
    #     return self._gas

    # @gas.setter
    # def gas(self, value):
    #     self._gas = SchDef.check_week_validity(value, 'gas')


        """
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
        q_elec: Number for maximum electrical plug process load [W/m2].
        q_gas: Number for maximum gas process load per unit area [W/m2].
        q_light: Number for maximum light process load per unit area [W/m2].
        n_occ: Number for maximum number of occupants per unit area [person/m2].
        vent: Number for maximum ventilation rate per unit area [m3/s/m2].
        v_swh: Number for maximum hot water rate per unit area [L/hr/m2].
        bldtype: Number between 0 and 15 corresponding to the following building
            types: FullServiceRestaurant (0), Hospital (1), LargeHotel (2),
            LargeOffice (3), MediumOffice (4), MidRiseApartment (5), OutPatient (6),
            PrimarySchool (7), QuickServiceRestaurant (8), SecondarySchool (9),
            SmallHotel (10), SmallOffice (11), StandaloneRetail (12), StripMall (13),
            SuperMarket (14), Warehouse (15). Additional building types can be defined
            with a number greater then 15. This value is used to reference the fraction
            of urban area the SchDef object defines in the UWG bld matrix.
        builtera: Number between 0 and 2 corresponding to the following built eras:
            Pre-1980s (0), Post1980s (1), New construction (2). This value is used to
            reference the fraction of urban area the SchDef object defines in the UWG
            bld matrix.
    """

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
            "elec": _example_week,
            "gas": _example_week,
            "light": _example_week,
            "occ": _example_week,
            "cool": _example_week,
            "heat": _example_week,
            "swh": _example_week,
            "q_elec" = q_elec,
            "q_gas" = q_gas,
            "q_light" = q_light,
            "n_occ" = n_occ,
            "vent" = vent,
            "v_swh" = v_swh,
            "bldtype": 0,
            "builtera": 1
            }
        """
        assert data['type'] == 'SchDef', 'Expected ' \
            'SchDef dictionary. Got {}.'.format(data['type'])

        return cls(elec=data['elec'], gas=data['gas'], light=data['light'],
                   occ=data['occ'], cool=data['cool'], heat=data['heat'],
                   swh=data['swh'], q_elec=data['q_elec'], q_gas=data['q_gas'],
                   q_light=data['q_light'], n_occ=data['n_occ'], vent=data['vent'],
                   v_swh=data['v_swh'], bldtype=data['bldtype'],
                   builtera=data['builtera'])

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
        base['q_elec'] = self.q_elec
        base['q_gas'] = self.q_gas
        base['q_light'] = self.q_light
        base['n_occ'] = self.n_occ
        base['vent'] = self.vent
        base['v_swh'] = self.v_swh
        base['bldtype'] = self.bldtype
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
                    'must contain 3 lists of numbers. Got : {}.'.format(name, val)
        return week

    def __repr__(self):
        return 'Schedule, bldtype: {}\n builtera: {}\n q_elec: {}\n q_gas: {}\n ' \
            'q_light: {}\n n_occ: {}\n vent: {}\n v_swh: {}\n elec: {}\n ' \
            'gas: {}\n light: {}\n occ: {}\n cool: {}\n heat: {}\n swh: {}\n'.format(
                self.bldtype, self.builtera, self.q_elec, self.q_gas, self.q_light,
                self.n_occ, self.vent, self.v_swh, self.elec, self.gas, self.light,
                self.occ, self.cool, self.heat, self.swh)
