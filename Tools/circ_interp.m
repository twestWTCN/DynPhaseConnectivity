function y2 = circ_interp(x,y1,xx)
y2 = [y1 y1 y1];
xx2 = unwrap([xx xx xx]);
x2  = unwrap([x x x]);
y2 = spline(x2,y2,xx2);
y2 = y2(length(xx)+1:(length(xx)+1)+length(xx)-1);