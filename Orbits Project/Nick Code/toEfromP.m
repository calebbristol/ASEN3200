function DCM = toEfromP(OMEGA,i,omega)
    cO = cos(OMEGA);
    sO = sin(OMEGA);
    ci = cos(i);
    si = sin(i);
    co = cos(omega);
    so = sin(omega);
    DCM = [co*cO-so*ci*sO co*sO+so*ci*cO so*si;...
            -so*cO-co*ci*sO -so*sO+co*ci*cO co*si;...
            si*sO -si*cO ci]';
end