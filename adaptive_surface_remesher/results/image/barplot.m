x = 10:10:60;
y{1} = [0.20       6.09      20.89      27.85      28.04      16.93];
y{2} = [0.32       1.81       4.40       9.92      31.38      52.17];
y{3} = [0.05       0.68       1.73       4.56      29.77      63.22];
y{4} = [0.04       0.32       0.89       2.92      27.99      67.84];
y{5} = [0.05       0.25       0.68       2.14      27.76      69.14];
y{6} = [0.02       0.16       0.57       2.19      26.99      70.07];

for i = 1:6
bar(x,y{i},"r"); axis([0 70 0 100])
hx=xlabel("Minimum Angle");
hy=ylabel("% of triangles");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);
print(["mp" sprintf("%d", i)  ".eps"],"-deps","-FArial", "-tight", "-color");
endfor
