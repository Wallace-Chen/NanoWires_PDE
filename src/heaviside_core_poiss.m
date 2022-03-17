% this function returns a value of 1 for inputted x and y values that lie
% within the inner core region of the nanowire
function inside = heaviside_core_poiss(core_side_length, x, y)

    global debug;
    
    global avg_side_length_poiss;
    global degree_of_polygon;
    
    % fudge factor so that all points along the boundary are accounted for
    delta = avg_side_length_poiss / 5;
    
    %% parameter check
    if degree_of_polygon ~= 3 && degree_of_polygon ~= 6 && degree_of_polygon ~= 0
       disp("Error: heaviside_core_poiss, the geometry not supported!");
       return;
    end
    
    %% check if (x,y) is within the region    
    if degree_of_polygon == 0
        inside = (x^2 + y^2 <= (core_side_length + delta)^2);
    elseif degree_of_polygon == 3
        inside = (y >= -core_side_length*sqrt(3)/6 - delta).*...
                 (y <= sqrt(3)*x+core_side_length*sqrt(3)/3 + 2*delta).*...
                 (y <= -sqrt(3)*x+core_side_length*sqrt(3)/3 + 2*delta);
    elseif degree_of_polygon == 6
        inside = (y <= sqrt(3)/2*core_side_length + delta).*...
                 (y <= -sqrt(3)*(x-core_side_length) + 2*delta).*...
                 (y >= sqrt(3)*(x-core_side_length) - 2*delta).*...
                 (y >= -sqrt(3)/2*core_side_length - delta).*...
                 (y >= -sqrt(3)*(x+core_side_length) - 2*delta).*...
                 (y <= sqrt(3)*(x+core_side_length) + 2*delta);
    end
    
    %% debug purpose
    if debug
       % too much infos flush the output
       % disp(fprintf("debug: heaviside_core_poiss, (%f, %f) is reported to be within the region with length %f, degree %d: %d\n", x, y, core_side_length, degree_of_polygon, inside)); 
    end
    
end