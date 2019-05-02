function be_t = be_t_calculate(ac_ratio, at_ratio, t)
a = at_ratio*t;
be_t = 1-(pi*a)/(2*t*(2+ac_ratio/at_ratio));
