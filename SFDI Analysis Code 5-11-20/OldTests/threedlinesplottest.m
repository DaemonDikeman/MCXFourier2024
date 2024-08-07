P0 = rand(7,3) ;
P1 = rand(7,3) ;
X = [P0(:,1) P1(:,1)] ;
Y = [P0(:,2) P1(:,2)] ;
Z = [P0(:,3) P1(:,3)] ;
plot3(X',Y',Z', '.')
hold on
% plot3(X',Y',Z','.')
hold off

%%
Wf=Xf.*Yf;