clear all; close all; clc;
A = load("test_output.dat");

x = A(:,1);
y = A(:,2);

for ii = 1:length(x)
    plot(x(ii),y(ii),'o'); hold on
    text(x(ii),y(ii),int2str(ii));
end