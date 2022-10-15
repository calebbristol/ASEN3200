function [num_launches, num_spacecraft, satellite_list] = loadConstellation(filename)

num_launches = 0;
num_spacecraft = 0;
satellite_list.name = '';
satellite_list.oe0 = NaN(6,1);

%1) extract the constellation structure from the json file
fid = fopen(filename);
raw = fread(fid,inf);
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

%2) read all of the launches and payloads to understand how many launches
% and spacecraft are in the constellation; note, this will be useful in
% Part 2!

num_launches = length(data.launches);
for i=1:num_launches
    num_spacecraft = num_spacecraft + length(data.launches(i).payload);
end

%3) RECOMMENDED: Pre-allocate the satellite_list struct
satellite_list(num_spacecraft).name = '';
satellite_list(num_spacecraft).oe0 = NaN(6,1);

%4) Populate each entry in the satellite struct list with its name and
%initial orbit elements [a,e,i,Om,om,f] at time t0
i = 1;
for j=1:num_launches
    oe0Temp1 = cell2mat(struct2cell((data.launches(j).orbit)));
    for k=1:length(data.launches(j).payload)
        nameTemp = data.launches(j).payload(k).name;
        oe0Temp2 = data.launches(j).payload(k).f;
        satellite_list(i).name = nameTemp;
        satellite_list(i).oe0 = [oe0Temp1; oe0Temp2];
        i = i + 1;
    end
end
