classdef RSMDef
    % Rural Site & Vertical Diffusion Model (VDM)
    % Calculates the vertical profiles of air temperature above the weather
    % station per 'The UWG' (2012) Eq. 4, 5, 6.
    
    properties
        lat;           % latitude (deg)
        lon;           % longitude (deg)
        GMT;           % GMT hour correction 
        height         % average obstacle height (m)
        z0r;           % rural roughness length (m)
        disp;          % rural displacement length (m)
        z;             % vertical height (m)
        dz;            % vertical discretization (m)
        nz0;           % layer number at zmt (m)
        nzref;         % layer number at zref (m)
        nzfor;         % layer number at zfor (m)
        nz10;          % layer number at zmu (m)
        nzi;           % layer number at zi_d (m)
        tempProf;      % potential temperature profile at the rural site (K)
        presProf;      % pressure profile at the rural site (Pa)
        tempRealProf;  % real temperature profile at the rural site (K)
        densityProfC;  % density profile at the center of layers (kg m-3)
        densityProfS;  % density profile at the sides of layers (kg m-3)
        windProf;      % wind profile at the rural site (m s-1)
        ublPres;       % Average pressure at UBL (Pa)
    end
    
    methods
        function obj = RSMDef(lat,lon,GMT,height,T_init,P_init,parameter)
            % class constructor
            load -ascii z_meso.txt;
            if(nargin > 0)
                obj.lat = lat;
                obj.lon = lon;
                obj.GMT = GMT;
                obj.height = height;
                obj.z0r = 0.1*height;
                obj.disp = 0.5*height;
                % vertical grid at the rural site
                obj.z  = zeros(numel(z_meso)-1,1);
                obj.dz = zeros(numel(z_meso)-1,1);
                for zi=1:numel(z_meso)-1
                    obj.z(zi) = 0.5*(z_meso(zi)+z_meso(zi+1));
                    obj.dz(zi) = z_meso(zi+1) - z_meso(zi);
                end
                ll = 1;
                mm = 1;
                nn = 1;
                oo = 1;
                pp = 1;
                for iz=1:55  
                   if (obj.z(iz)>=parameter.tempHeight && ll==1) 
                      obj.nz0 = iz;
                      ll = 0;
                   end  
                   if (obj.z(iz)>=parameter.refHeight && mm==1) 
                      obj.nzref = iz;
                      mm = 0;
                   end   
                   if (obj.z(iz)>=parameter.nightBLHeight && nn==1) 
                      obj.nzfor = iz;
                      nn = 0;
                   end   
                   if (obj.z(iz)>=parameter.windHeight && oo==1) 
                      obj.nz10 = iz;
                      oo = 0;
                   end   
                   if (obj.z(iz)>=parameter.dayBLHeight && pp==1) 
                      obj.nzi = iz;
                      pp = 0;
                   end   
                end
                % vertical profiles at the rural site
                obj.tempProf = ones(1,obj.nzref)*T_init;
                obj.presProf = ones(1,obj.nzref)*P_init;
                for iz=2:obj.nzref;
                   obj.presProf(iz) = (obj.presProf(iz-1)^(parameter.r/parameter.cp)-...
                       parameter.g/parameter.cp*(P_init^(parameter.r/parameter.cp))*(1./obj.tempProf(iz)+...
                       1./obj.tempProf(iz-1))*0.5*obj.dz(iz))^(1./(parameter.r/parameter.cp));
                end
                obj.tempRealProf = ones(1,obj.nzref)*T_init;
                for iz=1:obj.nzref;
                   obj.tempRealProf(iz)=obj.tempProf(iz)*...
                       (obj.presProf(iz)/P_init)^(parameter.r/parameter.cp);
                end
                obj.densityProfC = ones(1,obj.nzref);
                for iz=1:obj.nzref;
                   obj.densityProfC(iz)=obj.presProf(iz)/parameter.r/obj.tempRealProf(iz);
                end
                                obj.densityProfS = obj.densityProfC(1)*ones(1,obj.nzref+1);
                for iz=2:obj.nzref;
                   obj.densityProfS(iz)=(obj.densityProfC(iz)*obj.dz(iz-1)+...
                       obj.densityProfC(iz-1)*obj.dz(iz))/(obj.dz(iz-1)+obj.dz(iz));
                end
                obj.densityProfS(obj.nzref+1)=obj.densityProfC(obj.nzref);
                obj.windProf = ones(1,obj.nzref);
            end
        end
            
        % Ref: The UWG (2012), Eq. (4)
        function obj = VDM(obj,forc,rural,parameter,simTime)

            obj.tempProf(1) = forc.temp;    % Lower boundary condition  
            % compute pressure profile
            for iz=obj.nzref:-1:2
               obj.presProf(iz-1)=(obj.presProf(iz)^(parameter.r/parameter.cp)+...
                   parameter.g/parameter.cp*(forc.pres^(parameter.r/parameter.cp))*...
                   (1./obj.tempProf(iz)+1./obj.tempProf(iz-1))*...
                   0.5*obj.dz(iz))^(1./(parameter.r/parameter.cp));
            end
            % compute the real temperature profile
            for iz=1:obj.nzref
               obj.tempRealProf(iz)=obj.tempProf(iz)*...
                   (obj.presProf(iz)/forc.pres)^(parameter.r/parameter.cp);
            end
            % compute the density profile
            for iz=1:obj.nzref
               obj.densityProfC(iz)=obj.presProf(iz)/parameter.r/obj.tempRealProf(iz);
            end
            obj.densityProfS(1)=obj.densityProfC(1);
            for iz=2:obj.nzref
               obj.densityProfS(iz)=(obj.densityProfC(iz)*obj.dz(iz-1)+...
                   obj.densityProfC(iz-1)*obj.dz(iz))/(obj.dz(iz-1)+obj.dz(iz));
            end
            obj.densityProfS(obj.nzref+1)=obj.densityProfC(obj.nzref);
            
            % Ref: The UWG (2012), Eq. (5)
            % compute diffusion coefficient
            [cd,ustarRur] = DiffusionCoefficient(obj.densityProfC(1),...
                obj.z,obj.dz,obj.z0r,obj.disp,...
                obj.tempProf(1),rural.sens,obj.nzref,forc.wind,...
                obj.tempProf,parameter);
            % solve diffusion equation
            obj.tempProf = DiffusionEquation(obj.nzref,simTime.dt,...
                obj.tempProf,obj.densityProfC,obj.densityProfS,cd,obj.dz);
            % compute wind profile
            for iz=1:obj.nzref
                obj.windProf(iz) = ustarRur/parameter.vk*...
                    log((obj.z(iz)-obj.disp)/obj.z0r);
            end
            % Average pressure
            obj.ublPres = 0;
            for iz=1:obj.nzfor
                obj.ublPres = obj.ublPres +...
                    obj.presProf(iz)*obj.dz(iz)/...
                    (obj.z(obj.nzref)+obj.dz(obj.nzref)/2);
            end
        end
    end
end

function co = DiffusionEquation(nz,dt,co,da,daz,cd,dz)
    % Reference?

    cddz = zeros(nz+2,1);
    a = zeros(nz,3);
    c = zeros(nz,1);
    %--------------------------------------------------------------------------
    cddz(1)= daz(1)*cd(1)/dz(1);
    for iz=2:nz
       cddz(iz) = 2.*daz(iz)*cd(iz)/(dz(iz)+dz(iz-1));
    end
    cddz(nz+1) = daz(nz+1)*cd(nz+1)/dz(nz);
    %--------------------------------------------------------------------------
    a(1,1)=0.;
    a(1,2)=1.;
    a(1,3)=0.;
    c(1)=co(1);       
    for iz=2:nz-1
       dzv=dz(iz);
       a(iz,1)=-cddz(iz)*dt/dzv/da(iz);
       a(iz,2)=1+dt*(cddz(iz)+cddz(iz+1))/dzv/da(iz);
       a(iz,3)=-cddz(iz+1)*dt/dzv/da(iz);
       c(iz)  =co(iz);             
    end    
    a(nz,1)=-1.;
    a(nz,2)=1.;
    a(nz,3)=0.;
    c(nz)  =0;
    %--------------------------------------------------------------------------
    co = Invert (nz,a,c);
                                                   
end

function [Kt,ustar] = DiffusionCoefficient(rho,z,dz,z0,disp,...
    tempRur,heatRur,nz,uref,th,parameter)

    % Initialization
    Kt = zeros(1,nz+1);
    ws = zeros(1,nz);
    te = zeros(1,nz);
    % Friction velocity (Louis 1979)
    ustar = parameter.vk*uref/log((10.-disp)/z0);
    % Monin-Obukhov length
    lengthRur = max(- rho*parameter.cp*ustar^3*tempRur/parameter.vk/parameter.g/heatRur,-50);
    % Unstable conditions
    if gt(heatRur,1e-2)
        % Convective velocity scale
        wstar = (parameter.g*heatRur*parameter.dayBLHeight/rho/parameter.cp/tempRur)^(1/3);
        % Wind profile function
        phi_m = (1-8.*0.1*parameter.dayBLHeight/lengthRur)^(-1./3.);
        for iz=1:nz
            % Mixed-layer velocity scale
            ws(iz) = (ustar^3+phi_m*parameter.vk*wstar^3*z(iz)/parameter.dayBLHeight)^(1./3.);
            % TKE approximation
            te(iz) = max(ws(iz)^2.,0.01);
        end 
    % Stable and neutral conditions
    else
        for iz=1:nz
            % TKE approximation
            te(iz) = max(ustar^2.,0.01);
        end 
    end
    % lenght scales (l_up, l_down, l_k, l_eps)
    [dlu,dld] = DissipationBougeault(parameter.g,nz,z,dz,te,th);
    [dld,dls,dlk]= LengthBougeault(nz,dld,dlu,z);
    % Boundary-layer diffusion coefficient
    for iz=1:nz
       Kt(iz) = 0.4*dlk(iz)*sqrt(te(iz));
    end
    Kt(nz+1) = Kt(nz);

end

function [dlu,dld] = DissipationBougeault(g,nz,z,dz,te,pt)

    dlu = zeros(nz,1);
    dld = zeros(nz,1);
    for iz=1:nz
        zup=0.;
        dlu(iz)=z(nz+1)-z(iz)-dz(iz)/2.;
        zzz=0.;
        zup_inf=0.;
        beta=g/pt(iz);
        for izz=iz:nz-1
           dzt=(dz(izz+1)+dz(izz))/2.;
           zup=zup-beta*pt(iz)*dzt;
           zup=zup+beta*(pt(izz+1)+pt(izz))*dzt/2.;
           zzz=zzz+dzt;
           if (lt(te(iz),zup) && ge(te(iz),zup_inf))
             bbb=(pt(izz+1)-pt(izz))/dzt;
             if ne(bbb,0)
                tl=(-beta*(pt(izz)-pt(iz))+...
                sqrt( max(0.,(beta*(pt(izz)-pt(iz)))^2.+...
                2.*bbb*beta*(te(iz)-zup_inf))))/bbb/beta;
             else
             tl=(te(iz)-zup_inf)/(beta*(pt(izz)-pt(iz)));
             end            
             dlu(iz)=max(1.,zzz-dzt+tl);
           end
           zup_inf=zup;
        end
        zdo=0.;
        zdo_sup=0.;
        dld(iz)=z(iz)+dz(iz)/2.;
        zzz=0.;
        for izz=iz:-1:2
            dzt=(dz(izz-1)+dz(izz))/2.;
            zdo=zdo+beta*pt(iz)*dzt;
            zdo=zdo-beta*(pt(izz-1)+pt(izz))*dzt/2.;
            zzz=zzz+dzt;
            if (lt(te(iz),zdo) && ge(te(iz),zdo_sup))
               bbb=(pt(izz)-pt(izz-1))/dzt;
               if ne(bbb,0.)
                 tl=(beta*(pt(izz)-pt(iz))+...
                 sqrt( max(0.,(beta*(pt(izz)-pt(iz)))^2.+...
                 2.*bbb*beta*(te(iz)-zdo_sup))))/bbb/beta;
               else
                 tl=(te(iz)-zdo_sup)/(beta*(pt(izz)-pt(iz)));
               end
               dld(iz)=max(1.,zzz-dzt+tl);
            end
            zdo_sup=zdo;
        end
    end  
end

function [dld,dls,dlk] = LengthBougeault(nz,dld,dlu,z)

    dlg = zeros(nz,1);
    dls = zeros(nz,1);
    dlk = zeros(nz,1);
    for iz=1:nz
        dlg(iz)=(z(iz)+z(iz+1))/2.;
    end

    for iz=1:nz
        dld(iz)=min(dld(iz),dlg(iz));
        dls(iz)=sqrt(dlu(iz)*dld(iz));
        dlk(iz)=min(dlu(iz),dld(iz));         
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
