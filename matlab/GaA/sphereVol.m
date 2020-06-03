function sV = sphereVol(r,n)

sV = r.^n.*pi.^(n./2)./gamma(n./2+1);