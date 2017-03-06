function h=diva_hanning(N1,N2,N3)
h=[zeros(1,N1),.5+.5*cos(pi*(-N2:0)/(N2+1)),.5+.5*cos(pi*(1:N3)/(N3+1))];
h=h/sum(h);