size(4.2cm,2cm);

pair A=(0,0), B=(4.2,0), C=(4.2,2), D=(0,2);


draw (A--B--C--D--cycle);
//draw("$\partial_x C = 0$",(A+B)*0.5,W);
draw("Diagonal radius, $r_d$",A+(0.4,0.5),E);
draw("Axial radius, $r_h$",A+(0.4,1.5),E);
draw(box((0.2,0.4),(0.4,0.6)));

draw(shift(0.2,1.4)*scale(0.2)*(A--(A+(0.5,1))--(A+(1,0))--cycle));
//draw("$\partial_y C = 0$",(C+B)*0.5,N);
//draw("$\partial_y C = 0$",(A+D)*0.5,S);

