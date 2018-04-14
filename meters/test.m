clear all; close all; clc;
A = load("test_output.dat");

x = A(:,1);
y = A(:,2);
z = A(:,3);

for ii = 1:length(x)
    plot3(x(ii),y(ii),z(ii),'o'); hold on
    text(x(ii),y(ii),z(ii),int2str(ii));
end