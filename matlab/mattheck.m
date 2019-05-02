function M=mattheck(ac_ratio, at_ratio)

F1 = -1.9071*at_ratio + 1.5151*ac_ratio.^0.16596*at_ratio.^2-...
    21.52*ac_ratio.^2.1419*at_ratio.^3+0.34216*at_ratio.^4;
F2 = -0.74+3.855*ac_ratio-3.825*ac_ratio.^2-2.89*ac_ratio.^3+...
    4.356*ac_ratio.^4;
F3 = 1-at_ratio.^1.4;
M = (1-F1.*F2).*F3;
