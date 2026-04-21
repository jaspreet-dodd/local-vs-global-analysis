function y = softclip(x,k,r)
%
% x: data
% k: softclip hardness
% r: clip range

z = 2.*x./diff(r) - 1;

y = sign(z).*atan(abs(z).^k).^(1/k);
y = (y + 1).*diff(r)./2;