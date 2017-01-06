% N matrix
function n=Ne(xi, eta)
  n = 1/4*[ (1-xi)*(1-eta) 0 (1+xi)*(1-eta) 0 (1+xi)*(1+eta) 0 (1-xi)*(1+eta) 0; 
            0 (1-xi)*(1-eta) 0 (1+xi)*(1-eta) 0 (1+xi)*(1+eta) 0 (1-xi)*(1+eta); ]; 
end