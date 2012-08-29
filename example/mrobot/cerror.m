t1 = load("base-theta.true.out");
t2 = load("base-theta.out");
x1 = load("base-x.true.out");
x2 = load("base-x.out");
z1 = load("base-z.true.out");
z2 = load("base-z.out");

terr_p = sum((t1(10000:50000,1) - t2(10000:50000,1)).^2)/40000
terr_d = sum((t1(10000:50000,2) - t2(10000:50000,2)).^2)/40000
xerr_p = sum((x1(10000:50000,1) - x2(10000:50000,1)).^2)/40000
xerr_d = sum((x1(10000:50000,2) - x2(10000:50000,2)).^2)/40000
zerr_p = sum((z1(10000:50000,1) - z2(10000:50000,1)).^2)/40000
zerr_d = sum((z1(10000:50000,2) - z2(10000:50000,2)).^2)/40000

fprintf(1,"Mean positional error: %f\n", mean([terr_p xerr_p zerr_p]));
fprintf(1,"Mean derivative error: %f\n", mean([terr_d xerr_d zerr_d]));

