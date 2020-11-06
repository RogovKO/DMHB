function [omega_s, T_s, amp_s, fi_s, DC_s] = osciprof(dde_z,time,n)
%Determine Oscillatory profile from DATA
%   Detailed explanation goes here

int=round(length(time')/1.5);

len = length(dde_z(:,int:end)');
NFFT = 2^nextpow2(len);

FF = fft(dde_z(:,int:end)',NFFT)/len;

DC_s = FF(1,:);
DC_s = DC_s';

% amp_s = sqrt(mean(dde_z(:,int:end)'.^2));

[amp_s, loc] = max(dde_z(:,int:end)');

amp_s = amp_s'-DC_s;

[~, loc] = findpeaks(dde_z(1,:));
T_s = (time(loc(end))-time(loc(end-10)))/10;
omega_s = 2*pi/T_s;

it = 1;
fi_s=[];

while it <= n
    PhDiff = phdiffmeasure(dde_z(1,(int:end)), dde_z(it,(int:end)));
    fi_s = [fi_s;PhDiff];
    it = it + 1;
end

fir_s= rad2deg(fi_s);

di1 = ['Frequency is equal to ',num2str(omega_s)];
di2 = ['Period is equal to ',num2str(T_s)];


disp(di1)
disp(di2)


disp('DC is ');
disp(DC_s)

disp('Amplitudes are equal to ');
disp(amp_s)

disp('The phase shift is (RAD)');
disp(fi_s)
 
disp('The phase shift is (DEG)');
disp(fir_s)



end

