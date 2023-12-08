function f = integrand_conv(xvec, x, a1, b1, a2, b2)
    i1 = asymstab(xvec, a1, b1);
    i2 = asymstab(x-xvec, a2, b2);
    f = (i1.*i2)';
end