CL_EPW_PATH = 'C:\Sim\UWG4.1\data\';
CL_EPW = 'rural_weather_data_cap.epw';

CL_RE_PATH = 'C:\sim\UWG4.1\output';
CL_RE = 'rural_weather_data_cap_UWG.epw';

CL_XML_PATH = 'C:\sim\UWG4.1\data';
% CL_XML = {
%     'BackBayStation_27.xml'
%      };
% CL_XML = {
%     'initialize.m'
%      };
CL_XML = {
    'RunUWG_CAP.xlsm'
     };
 
for i = 1:length(CL_XML)
    [new_climate_file] = UWG(CL_EPW_PATH,CL_EPW,CL_XML_PATH,CL_XML{i},CL_RE_PATH,CL_RE);
end