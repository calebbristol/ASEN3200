function inLoS = testLoS(r_site,r_sc,elevation_limit)

foundAngle = dot(r_site,r_sc)/(norm(r_site)*norm(r_sc));

if foundAngle > elevation_limit
    inLoS = 1;
else
    inLoS = 0;
end