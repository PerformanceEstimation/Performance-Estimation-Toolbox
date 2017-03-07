function [ g ] = subgradient( x,f )
[g,~]=f.oracle(x);
end

