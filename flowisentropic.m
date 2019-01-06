function fval = flowisentropic (M)
  G = 1.4;
  A_ratio = 2.45;
  
  fval = ((1/M)*(2/(G+1) + ((G-1)/(G+1))*M^2)^(G+1)/((2*(G-1)))) - A_ratio;
  
endfunction
