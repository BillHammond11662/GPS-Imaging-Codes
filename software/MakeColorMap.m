
vumin = -3;
vumax = 3;

NC=64;
z1 = linspace(0,1,NC/2)';
z1f = flipud(z1);
cmap = [ones(NC/2,1) z1 z1; 
        z1f          z1f  ones(NC/2,1)];
cmap = flipud(cmap);
colormap(cmap);
m = (NC - 1)/(vumax-vumin);
b = 1 - m*vumin;
