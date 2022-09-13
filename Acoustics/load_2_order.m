function [l, beta, r, Lphi60] = load_2_order(dr, deltaphi)

fid = fopen('harmonics_load.txt');
nharmonics = fgetl(fid);        %read total number of harmonics
nharmonics = str2num(nharmonics);
nharmonics = linspace(1,nharmonics,nharmonics);
radiuses=split(fgetl(fid),';'); %read radius coordinate for each harmonic

%extract coefficients data
L = zeros(length(radiuses), length(nharmonics) + 1);
M = zeros(length(radiuses), length(nharmonics) );
fclose(fid);

harmonics = importdata('harmonics_load.txt','\t',3).data;
z = importdata('harmonics_load.txt','\t',3);
i = 0;
k = 0;
for i=1:length(radiuses)
    for j= 1:length(nharmonics)+1
        L(i,j) = importdata('harmonics_load.txt','\t',3).data(k+j,1);
        M(i,j) = importdata('harmonics_load.txt','\t',3).data(k+j,2);
    end
    k = k + length(nharmonics)+1;
end
M(:,1) = [];
L = L';
M = M';




fid = fopen('harmonics_beta.txt');
nharmonics = fgetl(fid)        %read total number of harmonics
nharmonics = str2num(nharmonics);
nharmonics = linspace(1,nharmonics,nharmonics);
radiuses=split(fgetl(fid),';'); %read radius coordinate for each harmonic

%extract coefficients data
L_beta = zeros(length(radiuses), length(nharmonics) + 1);
M_beta = zeros(length(radiuses), length(nharmonics) );
fclose(fid);

harmonics = importdata('harmonics_beta.txt','\t',3).data;
z = importdata('harmonics_beta.txt','\t',3);
i = 0;
k = 0;
for i=1:length(radiuses)
    for j= 1:length(nharmonics)+1
        L_beta(i,j) = importdata('harmonics_beta.txt','\t',3).data(k+j,1);
        M_beta(i,j) = importdata('harmonics_beta.txt','\t',3).data(k+j,2);
    end
    k = k + length(nharmonics)+1;
end
M_beta(:,1) = [];
L_beta = L_beta';
M_beta = M_beta';


%end parsing


%sample of calculating an harmonic for azimuths X
r = [0.25,0.4,0.55,0.75,0.85,0.90,0.95];
phi = linspace(0,2*pi,360/deltaphi);
lphi = zeros(length(r)+2 , 360/deltaphi);
beta = zeros(1, 360/deltaphi);


for m=1:length(r)
    lphi(m+1,:) = Fseriesval(L(:,m),M(:,m),phi); % columb corresponds to r/R = 0.25
end
beta(1,:) = Fseriesval(L_beta(:,1),M_beta(:,1),phi); % columb corresponds to r/R = 0.25
%lphi

r = [0.093,0.25,0.4,0.55,0.75,0.85,0.90,0.95,1];
Lphi60 = lphi(:,60);
for m=1:length(phi)
plot(r*28,lphi(:,m))
hold on
end

%%%%%%%%%%%%%%%

Pr1 = zeros(3,length(phi));
Pr2 = zeros(3,length(phi));
Pr3 = zeros(3,length(phi));
Pr4 = zeros(3,length(phi));


for m=1:length(phi)
    
        Pr1(:,m) = polyfit([0.093; 0.25; 0.4 ],lphi(1:3,m),2);
        Pr2(:,m) = polyfit([0.4; 0.55; 0.75],lphi(3:5,m),2);
        Pr3(:,m) = polyfit([0.75; 0.85; 0.9],lphi(5:7,m),2);
        Pr4(:,m) = polyfit([0.9; 0.95; 1],lphi(7:9,m),2);

end


l1 = zeros(1/dr,length(phi));
l2 = zeros(1/dr,length(phi));
l3 = zeros(1/dr,length(phi));
l4 = zeros(1/dr,length(phi));



for m=1:length(phi)
    l1(:,m) = polyval( Pr1(:,m), linspace(0.093,0.4, 1/dr));
    l2(:,m) = polyval( Pr2(:,m), linspace(0.4,0.75, 1/dr));
    l3(:,m) = polyval( Pr3(:,m), linspace(0.75,0.9, 1/dr));
    l4(:,m) = polyval( Pr4(:,m), linspace(0.9,1, 1/dr));

end

l=[l1; l2; l3; l4];


r1=linspace(0.093,0.4, 1/dr);
r2=linspace(0.4,0.75, 1/dr);
r3= linspace(0.75,0.9, 1/dr);
r4= linspace(0.9,1, 1/dr);
r=[r1, r2, r3, r4];


end


