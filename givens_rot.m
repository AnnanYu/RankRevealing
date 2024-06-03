function rot = givens_rot(ap,aq)
if ap == 0
  s = 1; c = 0;
else
  dist2 = sqrt(ap^2+aq^2);
  c = abs(ap)/dist2;
  s= -sign(ap)*aq/dist2;
end
rot = [c,s;-s,c];
