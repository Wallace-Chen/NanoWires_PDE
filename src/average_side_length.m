% a function to compute the average side lenght of the triangular mesh
function ave = average_side_length(p, t)

    h3 = sqrt( (p(1,t(1,:)) - p(1,t(2,:))).^2 + (p(2,t(1,:)) - p(2,t(2,:))).^2 );
    h1 = sqrt( (p(1,t(2,:)) - p(1,t(3,:))).^2 + (p(2,t(2,:)) - p(2,t(3,:))).^2 );
    h2 = sqrt( (p(1,t(3,:)) - p(1,t(1,:))).^2 + (p(2,t(3,:)) - p(2,t(1,:))).^2 );
    
    ave = mean( [h1 h2 h3] );

end