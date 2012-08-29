% plots the controller performance for a given wheel
% example: plot_results("left-wheel-joint",1,10000)
%
% arguments:
% joint                the name of the joint to plot
% st                   the iteration to start plotting (optional)
% ed                   the iteration to stop plotting (optional)
%
% The first plot is the true and desired joint position.
% The second plot is the true and desired joint velocity.
% The third plot is the true and desired joint acceleration.
% The fourth plot is the set of motor torques, both composite and individual
%    (i.e., P, I, D, computed torque).

function plot_results(joint, st, ed)

file1 = [joint ".true.out"];
file2 = [joint ".out"];
file3 = [joint ".torques.out"];

% try to open the files for reading
[fid1,msg] = fopen(file1, "r");
[fid2,msg] = fopen(file2, "r");
[fid3,msg] = fopen(file3, "r");
if (fid1 != -1)
	fclose(fid1);
	x = load(file1);
else
	x = [];
end
if (fid2 != -1)
	fclose(fid2);
	y = load(file2);
else
	y = [];
end
if (fid3 != -1)
	fclose(fid3);
	z = load(file3);
else
	z = [];
end

% set st if not already set
if (nargin() == 1)
	st = 1
end

% set ed if not already set
if (nargin() < 3)
	if (size(x,1) != 0 && size(y,1) != 0)
		ed = min(size(x,1), size(y,1));
	elseif (size(x,1) != 0)
		ed = size(x,1);
	else
		ed = size(y,1);
	end
end

% plot the position (desired and true)
figure;
mx = inf;
mn = -inf;
if (size(x,1) > 0)
	mx_x = max(x(st:ed,1));
	mn_x = min(x(st:ed,1));
	if (size(y, 1) == 0)
		mx_y = inf;
		mn_y = -inf;
	else
		mx_y = max(y(st:ed,1));
		mn_y = min(y(st:ed,1));
	end
	mx = max(mx_x, mx_y);
	mn = min(mn_x, mn_y);
elseif (size(y,1) > 0)
	mx = max(y(st:ed,1));
	mn = min(y(st:ed,1));
end
if (mx - mn > 1e-4)
	axis([st/10000 ed/10000 max(-7,mn) min(7,mx)]);
else
	axis([st/10000 ed/10000 max(-7,mn)-1e-4 min(7,mx)+1e-4]);
end
hold;
if (size(x,1) > 0)
	plot(st/10000:.0001:ed/10000,x(st:ed,1), "r+;q (true);");
end
if (size(y,1) > 0)
	plot(st/10000:.0001:ed/10000,y(st:ed,1), "g*;q (des);");
end
xlabel("time");
hold off;

% plot the velocity (desired and true)
figure;
mx = inf;
mn = -inf;
if (size(x,1) > 0)
	mx_x = max(x(st:ed,2));
	mn_x = min(x(st:ed,2));
	if (size(y, 1) == 0)
		mx_y = inf;
		mn_y = -inf;
	else
		mx_y = max(y(st:ed,2));
		mn_y = min(y(st:ed,2));
	end
	mx = max(mx_x, mx_y);
	mn = min(mn_x, mn_y);
elseif (size(y,1) > 0)
	mx = max(y(st:ed,2));
	mn = min(y(st:ed,2));
end
if (mx - mn > 1e-4)
	axis([st/10000 ed/10000 max(-100,mn) min(100,mx)]);
else
	axis([ed/10000 ed/10000 max(-100,mn)-1e-4 min(100,mx)+1e-4]);
end
hold;
if (size(x,1) > 0)
	plot(st/10000:.0001:ed/10000,x(st:ed,2), "r+;qd (true);");
end
if (size(y,1) > 0)
	plot(st/10000:.0001:ed/10000,y(st:ed,2), "g*;qd (des);");
end
xlabel("time");
hold;

% plot the acceleration (desired and true)
figure;
mx = inf;
mn = -inf;
if (size(x,1) > 0)
	mx_x = max(x(st:ed,3));
	mn_x = min(x(st:ed,3));
	if (size(y, 1) == 0)
		mx_y = inf;
		mn_y = -inf;
	else
		mx_y = max(y(st:ed,3));
		mn_y = min(y(st:ed,3));
	end
	mx = max(mx_x, mx_y);
	mn = min(mn_x, mn_y);
elseif (size(y,1) > 0)
	mx = max(y(st:ed,3));
	mn = min(y(st:ed,3));
end

axis([st/10000 ed/10000 max(-100,mn) min(100,mx)]);
hold;
if (size(x,1) > 0)
	plot(st/10000:.0001:ed/10000, x(st:ed,3), "+;qdd (true);");
end
if (size(y,1) > 0)
	plot(st/10000:.0001:ed/10000, y(st:ed,3), ";qdd (des);");
end
xlabel("time");
hold;

% plot the torques
figure;
if (size(z,1) > 0)
	composite = sum(z(st:ed,1:4),2);
	mx = min(500, max(composite));
	mn = max(-500, min(composite));	
	axis([0 (ed-st+1)/10000 mn mx]);
	plot(st/10000:.0001:ed/10000,z(st:ed,1), ";invdyn;");
	hold;
	plot(st/10000:.0001:ed/10000,z(st:ed,2), ";P;");
	plot(st/10000:.0001:ed/10000,z(st:ed,3), ";D;");
	plot(st/10000:.0001:ed/10000,z(st:ed,4), ";I;");
	plot(st/10000:.0001:ed/10000,composite, ";composite;");
	hold;
end


