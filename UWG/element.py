
class Element(object):
    """
    UWG Element

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
    "-----------------------------------------\n"   +\
    "ERROR: the number of layer thickness must\n"   +\
    "match the number of layer materials\n"         +\
    "-----------------------------------------"
    def __init__(self, alb, emis, thicknessLst, materialLst, vegCoverage, T_init, horizontal,name=None):
        if len(thicknessLst) != len(materialLst):
            raise Exception(self.THICKNESSLST_EQ_MATERIALLST_MSG)
        else:
            self._name = name # purely for internal process
            self.albedo = alb
            self.emissivity = emis
            self.layerThickness = thicknessLst
            self.layerThermalCond = map(lambda z: 0, materialLst)
            self.layerVolHeat = map(lambda z: 0, materialLst)

            """
            for i = 1:numel(materialLst)
              self.layerThermalCond(i) = materialLst(i).thermalCond
              self.layerVolHeat(i) = materialLst(i).volHeat
            end
            """
            self.vegCoverage = vegCoverage
            """
            self.layerTemp = T_init*ones(numel(ThicknessLst),1)
            self.waterStorage = 0.
            self.infra = 0
            self.horizontal = horizontal
            self.sens = 0.
            """
    def __repr__(self):
        #returns some representative wall properties
        return "Element: {n}, depth={z}, e={a}, k={b}, Cp*dens={c}".format(
            n=self._name,
            z=sum(self.layerThickness),
            a=self.emissivity,
            b=str(sum(self.layerThermalCond)),
            c=str(sum(self.layerVolHeat))
            )

        """

        function obj = SurfFlux(obj,forc,parameter,simTime,humRef,tempRef,windRef,boundCond,intFlux)

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

        function t = Conduction(obj,dt,flx1,bc,temp2,flx2)
            t = obj.layerTemp;
            hc = obj.layerVolHeat;
            tc = obj.layerThermalCond;
            d = obj.layerThickness;
            % flx1  : net heat flux on surface
            % bc    : boundary condition parameter (1 or 2)
            % temp2 : deep soil temperature (ave of air temperature)
            % flx2  : surface flux (sum of absorbed, emitted, etc.)

            fimp=0.5;           % implicit coefficient
            fexp=0.5;           % explicit coefficient
            num = size(t,1);    % number of layers

            % mean thermal conductivity over distance between 2 layers
            tcp = zeros(num,1);
            % thermal capacity times layer depth
            hcp = zeros(num,1);
            % lower, main, and upper diagonals
            za = zeros(num,3);
            % RHS
            zy = zeros(num,1);
            %--------------------------------------------------------------------------
            hcp(1) = hc(1)* d(1);
            for j=2:num;
              tcp(j) = 2./(d(j-1)/tc(j-1)+d(j)/tc(j));
              hcp(j) = hc(j)*d(j);
            end
            %--------------------------------------------------------------------------
            za(1,1) = 0.;
            za(1,2) = hcp(1)/dt + fimp*tcp(2);
            za(1,3) = -fimp*tcp(2);
            zy(1) = hcp(1)/dt*t(1) - fexp*tcp(2)*(t(1)-t(2)) + flx1;
            %--------------------------------------------------------------------------
            for j=2:num-1;
              za(j,1) = fimp*(-tcp(j));
              za(j,2) = hcp(j)/dt+ fimp*(tcp(j)+tcp(j+1));
              za(j,3) = fimp*(-tcp(j+1));
              zy(j) = hcp(j)/dt*t(j)+fexp*(tcp(j)*t(j-1)-...
                  tcp(j)*t(j)-tcp(j+1)*t(j)+ tcp(j+1)*t(j+1));
            end
            %--------------------------------------------------------------------------
            if eq(bc,1) % het flux
                za(num,1) = fimp*(- tcp(num) );
                za(num,2) = hcp(num)/dt+ fimp* tcp(num);
                za(num,3) = 0.;
                zy(num) = hcp(num)/dt*t(num) + fexp*tcp(num)*(t(num-1)-t(num)) + flx2;
            elseif eq(bc,2) % deep-temperature
                za(num,1) = 0;
                za(num,2) = 1;
                za(num,3) = 0.;
                zy(num) = temp2;
            else
                disp('ERROR: check input parameters in the Conduction routine')
            end
            %--------------------------------------------------------------------------
            % zx=tridiag_ground(za,zb,zc,zy);
            zx = Invert(num,za,zy);
            t(:) = zx(:);

        end
     end
end

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

function x = Invert(nz,a,c)

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

end

"""
