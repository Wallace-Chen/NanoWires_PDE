% griddata_NaN has the same functionality as the griddata function already
% in MATLAB but ensures that it does not return any NaN numbers (which can
% cause errors for other programs)
function ZI = griddata_NaN(x, y, z, XI, YI)

    ZI = griddata(x, y, z, XI, YI, 'cubic');

    ZI( isnan(ZI) )=0;
    
end