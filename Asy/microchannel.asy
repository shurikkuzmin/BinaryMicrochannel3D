import three;
import graph3;
//currentprojection=orthographic(5,4,2,center=true);

currentprojection=perspective(
camera=(25.0851928432063,-30.3337528952473,19.3728775115443),
up=Z,
zoom=1,
autoadjust=false);

real r(real Theta, real Phi){return 1+0.05*(cos(Theta)*sin(Phi));}
triple f(pair z) {return r(z.x,z.y)*expi(z.x,z.y);}

size(5cm);
size3(5cm,5cm,5cm);

path3 face1=(-2,-1,0)--(2,-1,0)--(2,1,0)--(-2,1,0)--cycle;
path3 face2=(-2,-1,0)--(-2,-1,1)--(2,-1,1)--(2,-1,0)--cycle;
path3 face3=(-2,1,0)--(-2,1,1)--(2,1,1)--(2,1,0)--cycle;
path3 face4=(-2,-1,1)--(2,-1,1)--(2,1,1)--(-2,1,1)--cycle;

(Label(XY()*"$x$",align=-3Y));

//surface s=surface(f,(0,0),(pi,2pi),100,Spline);
draw("$z$",project(O--X),Arrow);
draw("$y$",project(O--Y),Arrow);
draw("$x$",project(O--Z),Arrow);
//draw(YZ);

label(project("W",Y,Z,(2,0,-0.3)));
label(project("L",X,Z,(0,-1,-0.3)));
label(project("H",Y,Z,(-2,-0.8,0.5)));


draw(face1);
draw(face2);
draw(face3);
draw(face4);
//draw(s);
