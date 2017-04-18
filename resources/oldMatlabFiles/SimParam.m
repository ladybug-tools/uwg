classdef SimParam
    %SIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt;            % Simulation time-step
        timeForcing;   % Weather data time-step
        month;         % Begin month
        day;           % Begin day of the month
        days;          % Number of days of simulation
        timePrint;     % time-step for printing outputs
        timeDay;       % number of timesteps in a design-day
        timeSim;       % number of timesteps of the simulation
        timeMax;       % total simulation time (s)
        nt;            % total number of timesteps
        timeFinal;     % final timestep of simulation
        timeInitial;   % initial timestep of simulation
        secDay;        % seconds of one day (s)
        hourDay;       % hour of the day (0 - 23hr)
        inobis;        % julian day at the end of each month
        julian;        % julian day
    end
    
    methods
        function obj = SimParam(dt,timefor,M,DAY,days)
            % class constructor
            if(nargin > 0)
                obj.dt = dt;
                obj.timeForcing = timefor;
                obj.month = M;
                obj.day = DAY;
                obj.days = days;
                obj.timePrint = timefor;
                obj.timeDay = 24*3600/timefor;
                obj.timeSim = obj.timeDay*days;
                obj.timeMax = 24.*3600.*days;
                obj.nt = fix(obj.timeMax/obj.dt+1);
                obj.inobis = [0,31,59,90,120,151,181,212,243,273,304,334];
                obj.julian = obj.inobis(M)+DAY-1;
                H1 = (obj.inobis(M)+DAY-1)*obj.timeDay;
                obj.timeInitial = H1 + 8;
                obj.timeFinal = H1 + obj.timeDay*obj.days - 1 + 8;
                obj.secDay = 0;
                obj.hourDay = 0;
            end
        end
        
        % update date
        function obj = UpdateDate(obj)
            obj.secDay = obj.secDay + obj.dt;
            if obj.secDay == 3600*24
                obj.day = obj.day+1;
                obj.julian = obj.julian+1;
                obj.secDay = 0;
                for j = 1:12
                    if obj.julian == obj.inobis(j)
                        obj.month = obj.month + 1;
                        obj.day = 1;
                    end
                end
            end
            obj.hourDay = floor(obj.secDay/3600);       % 0 - 23hr
        end
    end
end

