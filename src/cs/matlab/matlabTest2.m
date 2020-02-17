function vec_out = matlabTest(vec_in)

vec = detrend(vec_in);

[p,f] = pyulear(vec, 400, 2048, 100);

p = [p(2:end); p(2:end)];

xcor = xcorr(vec,p);

xcorhist = hist(xcor,10);

vec_out = resample(xcorhist,2048,10);

vec_out = vec_out(:);


end