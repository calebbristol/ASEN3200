function [num_launches, num_spacecraft, satellite_list] = loadConstellation(filename)
%DESCRIPTOIN: Ingests constellation description .json file and parses it
%into a list of structs with full initial orbit elements (km, s, rad) and
%satellite name.
%
%INPUTS:
% filename      A string indicating the name of the .json file to be parsed
%
%OUTPUTS:
% nl            Number of total launches
% ns            Total number of spacecraft between all launches
% satlist       Array of structs with 'name' and 'oe0' properties


%Temporary - just so the function runs the first time you use it.
%You'll need to change all of these!
num_launches = 0;
num_spacecraft = 0;
satellite_list.name = '';
satellite_list.oe0 = NaN(6,1);

%1) extract the constellation structure from the json file

fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);

%2) read all of the launches and payloads to understand how many launches
% and spacecraft are in the constellation; note, this will be useful in
% Part 2!
launches = val.launches;
num_launches = length(launches);
num_spacecraft = 0;
for i = 1:num_launches
    num_spacecraft = num_spacecraft + length(launches(i).payload);
end


%3) RECOMMENDED: Pre-allocate the satellite_list struct
satellite_list(num_spacecraft).name = '';
satellite_list(num_spacecraft).oe0 = NaN(6,1);

%4) Populate each entry in the satellite struct list with its name and
%initial orbit elements [a,e,i,Om,om,f] at time t0
icraft = 1;
for i = 1:num_launches
    for j = 1:length(launches(i).payload)
        satellite_list(icraft).name = launches(i).payload(j).name;
        satellite_list(icraft).oe0 = [launches(i).orbit.a ...
            launches(i).orbit.e launches(i).orbit.i ... 
            launches(i).orbit.Om launches(i).orbit.om ...
            launches(i).payload(j).f];
        icraft = icraft + 1;
    end
end


