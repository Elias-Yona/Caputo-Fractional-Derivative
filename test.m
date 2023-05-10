a = [1 1/16 4/5 3/2 1/25 6/5]; na = [3 2.5 2 1 0.5 0];
b = 1; nb = 0; t = [0:0.1:30]; u = 172/125*cos(4*t/5);
y0 = [1 4/5 -16/25]; y = sqrt(2)*sin(4*t/5+pi/4);
y1 = fode_caputo9(a,na,b,nb,y0,u,t,5);
max(abs(y-y1)), plot(t,y,t,y1)