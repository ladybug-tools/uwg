
class Element(object):
    """
    UWG Element

    # Note: In matlab not all instance variables are instantiated. They are assumed to be a 0-by-0 empty matrix
    # https://www.mathworks.com/help/matlab/matlab_oop/specifying-properties.html

    Attributes:
        albedo;          % outer surface albedo
        emissivity;      % outer surface emissivity
        layerThickness;  % vector of layer thicknesses (m)
        layerThermalCond;% vector of layer thermal conductivities (W m-1 K-1)
        layerVolHeat;    % vector of layer volumetric heat (J m-3 K-1)
        vegCoverage;     % surface vegetation coverage
        layerTemp;       % vector of layer temperatures (K)
        waterStorage;    % thickness of water film (m) (only for horizontal surfaces)
        horizontal;      % 1-horizontal, 0-vertical
        solRec;          % solar radiation received (W m-2)
        infra;           % net longwave radiation (W m-2)
        lat;             % surface latent heat flux (W m-2)
        sens;            % surface sensible heat flux (W m-2)
        solAbs;          % solar radiation absorbed (W m-2)
        aeroCond;        % convective heat transfer
        T_ext;           % external surface temperature
        T_int;           % internal surface temperature
        flux;            % external surface heat flux
    """

    THICKNESSLST_EQ_MATERIALLST_MSG = \
    "-----------------------------------------\n" +\
    "ERROR: the number of layer thickness must\n" +\
    "match the number of layer materials\n"
    "-----------------------------------------"
    CONDUCTION_INPUT_MSG = 'ERROR: check input parameters in the Conduction routine'

    def __init__(self, alb, emis, thicknessLst, materialLst, vegCoverage, T_init, horizontal,name=None):
        if len(thicknessLst) != len(materialLst):
            raise Exception(self.THICKNESSLST_EQ_MATERIALLST_MSG)
        else:
            self._name = name                                       # purely for internal process
            self.albedo = alb                                       # outer surface albedo
            self.emissivity = emis                                  # outer surface emissivity
            self.layerThickness = thicknessLst                      # vector of layer thicnesses (m)
            self.layerThermalCond = map(lambda z: 0, materialLst)   # vector of layer thermal conductivity (W m-1 K-1)
            self.layerVolHeat = map(lambda z: 0, materialLst)       # vector of layer volumetric heat (J m-3 K-1)

            #create list of layer k and (Cp*density) from materialLst properties
            for i in xrange(len(materialLst)):
              self.layerThermalCond[i] = materialLst[i].thermalCond
              self.layerVolHeat[i] = materialLst[i].volHeat

            self.vegCoverage = vegCoverage                          # surface vegetation coverage
            self.layerTemp = [T_init] * len(thicknessLst)           # vector of layer temperatures (K)
            self.waterStorage = 0.                                  # thickness of water film (m) for horizontal surfaces only
            self.infra = 0.                                         # net longwave radiation (W m-2)
            self.horizontal = horizontal                            # 1-horizontal, 0-vertical
            self.sens = 0.                                          # surface sensible heat flux (W m-2)

            # B/c we have to explicity define this in python. Set as None
            self.solRec = None                                       # solar radiation received (W m-2)
            self.lat = None                                          # surface latent heat flux (W m-2)
            self.solAbs = None                                       # solar radiation absorbed (W m-2)
            self.aeroCond = None                                     # convective heat transfer
            self.T_ext = None                                        # external surface temperature
            self.T_int = None                                        # internal surface temperature
            self.flux = None                                         # external surface heat flux

    def __repr__(self):
        #returns some representative wall properties
        s1 = "Element: {a}\n\tlayerNum={b}, totaldepth={c}\n\t".format(
            a=self._name,
            b=len(self.layerThickness),
            c=sum(self.layerThickness)
            )
        s2 = "e={d}, a={e}\n\tk_avg={f}, Cp*dens_avg={g}\n\tlayerTemp: {h}".format(
            d=self.emissivity,
            e=self.albedo,
            f=round(sum(self.layerThermalCond),2),
            g=round(sum(self.layerVolHeat),2),
            h=self.layerTemp
            )
        return s1 + s2

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def SurfFlux(self,forc,parameter,simTime,humRef,tempRef,windRef,boundCond,intFlux):
        pass
        """
        % Calculated per unit area (m^2)
        dens = forc.pres/(1000*0.287042*tempRef*(1.+1.607858*humRef)); % air density
        obj.aeroCond = 5.8+3.7*windRef;         % Convection coef (ref: UWG, eq. 12))

        if (obj.horizontal)     % For roof, mass, road

            % Evaporation (m s-1), Film water & soil latent heat
            if obj.waterStorage > 0
                qtsat = qsat(obj.layerTemp(1),forc.pres,parameter);
                eg = obj.aeroCond*parameter.colburn*dens*(qtsat-humRef)/parameter.waterDens/parameter.cp;
                obj.waterStorage = min(obj.waterStorage + simTime.dt*(forc.prec-eg),parameter.wgmax);
                obj.waterStorage = max(obj.waterStorage,0);
            else
                eg = 0;
            end
            soilLat = eg*parameter.waterDens*parameter.lv;

            % Winter, no veg
            if simTime.month < parameter.vegStart && simTime.month > parameter.vegEnd
                obj.solAbs = (1-obj.albedo)*obj.solRec;
                vegLat = 0;
                vegSens = 0;
            else    % Summer, veg
                obj.solAbs = ((1-obj.vegCoverage)*(1-obj.albedo)+...
                    obj.vegCoverage*(1-parameter.vegAlbedo))*obj.solRec;
                vegLat = obj.vegCoverage*parameter.grassFLat*(1-parameter.vegAlbedo)*obj.solRec;
                vegSens = obj.vegCoverage*(1.-parameter.grassFLat)*(1-parameter.vegAlbedo)*obj.solRec;
            end
            obj.lat = soilLat + vegLat;

            % Sensible & net heat flux
            obj.sens = vegSens + obj.aeroCond*(obj.layerTemp(1)-tempRef);
            obj.flux = - obj.sens+obj.solAbs+obj.infra-obj.lat;

        else     % Vertical surface (wall)
            obj.solAbs = (1-obj.albedo)*obj.solRec;
            obj.lat = 0;

            % Sensible & net heat flux
            obj.sens = obj.aeroCond*(obj.layerTemp(1)-tempRef);
            obj.flux = - obj.sens+obj.solAbs+obj.infra-obj.lat;
        end

        obj.layerTemp = Conduction(obj,simTime.dt,obj.flux,boundCond,forc.deepTemp,intFlux);
        obj.T_ext = obj.layerTemp(1);
        obj.T_int = obj.layerTemp(end);
    end
    """

    def Conduction(self, dt, flx1, bc, temp2, flx2):
        """
        Solve the conductance of heat based on of the element layers.
        arg:
            flx1  : net heat flux on surface
            bc    : boundary condition parameter (1 or 2)
            temp2 : deep soil temperature (ave of air temperature)
            flx2  : surface flux (sum of absorbed, emitted, etc.)

        key prop:
            za = [[ x00, x01, x02 ... x0w ]
                  [ x10, x11, x12 ... x1w ]
                            ...
                  [ xh0, xh1, xh2 ... xhw ]]

            where h = matrix row index    = element layer number
                  w = matrix column index = 3

        """
        t = self.layerTemp          # vector of layer temperatures (K)
        hc = self.layerVolHeat      # vector of layer volumetric heat (J m-3 K-1)
        tc = self.layerThermalCond  # vector of layer thermal conductivities (W m-1 K-1)
        d = self.layerThickness     # vector of layer thicknesses (m)

        # flx1                      : net heat flux on surface
        # bc                        : boundary condition parameter (1 or 2)
        # temp2                     : deep soil temperature (avg of air temperature)
        # flx2                      : surface flux (sum of absorbed, emitted, etc.)

        fimp = 0.5                  # implicit coefficient
        fexp = 0.5                  # explicit coefficient
        num = len(t)                # number of layers

        # Mean thermal conductivity over distance between 2 layers (W/mK)
        tcp = map(lambda tcon: 0, range(num))
        # Thermal capacity times layer depth (J/m2K)
        hcp = map(lambda tcap: 0, range(num))
        # lower, main, and upper diagonals
        za = map(lambda y_: map(lambda x_: 0, range(3)), range(num))
        # RHS
        zy = map(lambda rhs_: 0, range(num))

        #--------------------------------------------------------------------------
        # Define the column vectors for heat capactiy and conductivity
        hcp[0] = hc[0] * d[0]       # (J/m2K) = First row, define thermal capacity (J/m3K) * thickness (m)
        for j in xrange(1,num):     # From second row, define thermal conductivity and thermal capacity
            tcp[j] = 2. / (d[j-1] / tc[j-1] + d[j] / tc[j]) # (W/m2K) Mean of conductance (W/m2K) for layer j and j-1
            hcp[j] = hc[j] * d[j]   # (J/m2K)

        #--------------------------------------------------------------------------
        # Define the first row of za matrix, and RHS column vector
        za[0][0] = 0.
        za[0][1] = hcp[0]/dt + fimp*tcp[1]  # (J/m2K)/t + W/m2K
        za[0][2] = -fimp*tcp[1]             # W/m2K
        zy[0] = hcp[0]/dt*t[0] - fexp*tcp[1]*(t[0]-t[1]) + flx1 # W/m2 = K*(J/m2K)/t - K*W/m2K

        #--------------------------------------------------------------------------
        # ??? Define other rows
        for j in xrange(1,num-1):
          za[j][0] = fimp*(-tcp[j]) # W/m2K
          za[j][1] = hcp[j]/dt + fimp*(tcp[j]+tcp[j+1]) # J/m2K/t + W/m2K + W/m2K
          za[j][2] = fimp*(-tcp[j+1])   # W/m2K
          zy[j] = hcp[j]/dt * t[j] + fexp * \
            (tcp[j]*t[j-1] - tcp[j]*t[j] - tcp[j+1]*t[j] + tcp[j+1]*t[j+1]) # W/m2 = K*J/m2K/t * K*W/m2K - K*W/m2K - K*W/m2K + K*W/m2K

        #--------------------------------------------------------------------------
        # ??? Boundary conditions
        if self.is_near_zero(bc-1.): # heat flux
            za[num-1][0] = fimp * (-tcp[num-1])
            za[num-1][1] = hcp[num-1]/dt + fimp*tcp[num-1]
            za[num-1][2] = 0.
            zy[num-1] = hcp[num-1]/dt*t[num-1] + fexp*tcp[num-1]*(t[num-2]-t[num-1]) + flx2
        elif self.is_near_zero(bc-2.): # deep-temperature
            za[num-1][0] = 0.
            za[num-1][1] = 1.
            za[num-1][2] = 0.
            zy[num-1] = temp2
        else:
            raise Exception(self.CONDUCTION_INPUT_MSG)

        #--------------------------------------------------------------------------
        # zx=tridiag_ground(za,zb,zc,zy);
        zx = invert(num,za,zy)
        #t[:] = zx[:]
        #print '----fin----'

"""

function qsat = qsat(temp,pres,parameter)

    gamw = (parameter.cl - parameter.cpv) / parameter.rv;
    betaw = (parameter.lvtt/parameter.rv) + (gamw * parameter.tt);
    alpw = log(parameter.estt) + (betaw /parameter.tt) + (gamw *log(parameter.tt));
    work2 = parameter.r/parameter.rv;
    foes = zeros(size(temp));
    work1= zeros(size(temp));
    qsat = zeros(size(temp));
    for i=1:size(temp)
      % saturation vapor pressure
      foes(i) = exp( alpw - betaw/temp(i) - gamw*log(temp(i))  );
      work1(i)    = foes(i)/pres(i);
      % saturation humidity
      qsat(i) = work2*work1(i) / (1.+(work2-1.)*work1(i));
    end

end
"""

def invert(nz,a,c):
    """
    %--------------------------------------------------------------------------
    % Inversion and resolution of a tridiagonal matrix
    %          A X = C
    % Input:
    %  a(*,1) lower diagonal (Ai,i-1)
    %  a(*,2) principal diagonal (Ai,i)
    %  a(*,3) upper diagonal (Ai,i+1)
    %  c
    % Output
    %  x     results
    %--------------------------------------------------------------------------

    x = zeros(nz,1);

    for in=nz-1:-1:1
        c(in)=c(in)-a(in,3)*c(in+1)/a(in+1,2);
        a(in,2)=a(in,2)-a(in,3)*a(in+1,1)/a(in+1,2);
    end

    for in=2:nz
        c(in)=c(in)-a(in,1)*c(in-1)/a(in-1,2);
    end

    for in=1:nz
        x(in)=c(in)/a(in,2);
    end
    """
    return None#x
