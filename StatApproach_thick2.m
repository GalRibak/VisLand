Term = [151
 110
 65
 62
 74
 45
 49
 55
 82
 65
 66
 51
 46
 151
 90
 95
 104
 85
 126
 84
 97
 19
 87
 72
 264
 127
 79
 79
 116
];

lgt = [115
204
177
141
338
210
343
176
200
265
213
283
206
290
156
123
368
55
110
419
225
77
79
517
254
290
298
4
239
];

fpf =15; %flapping frequency in Hz
for i=1:29
   clearvars -except i Term lgt Dis MV0 MVt VT Res RES fpf Ps Es
    file = strcat('acth',num2str(i), '.csv');
    A = xlsread(file);
    
    [r,c] = size(A);

f=1000; 
    f_cutoff = 90; 
    fnorm = f_cutoff/(f/2);
    [b1,a1] = butter(4,fnorm,'low');
    A(:,1:3) = filtfilt(b1,a1,A(:,1:3));    
    A(:,4:6) = filtfilt(b1,a1,A(:,4:6));
    A(:,7:9) = filtfilt(b1,a1,A(:,7:9));

    A(:,10:12) = [A(:,1)-A(:,7),A(:,2)-A(:,8),A(:,3)-A(:,9)];
    A(:,13:15) = velocity(A(:,10:12),1000);
    A(:,16) = sqrt(A(:,13).^2 + A(:,14).^2);
    A(:,17) = sqrt(A(:,10).^2 + A(:,11).^2);
    Dis(i,1) = A(1,17); 
    fl(:,1) = and(A(:,17)>0.03, A(:,17)<0.04);
    fl(:,2) = and(A(:,17)>0.01, A(:,17)<0.02); % here 1-2
    A(:,18:20) = velocity([A(:,17) A(:,17) A(:,17)], 1000);
    A(:,16) = -A(:,18);
    V0=0;j=1; k=0;
    while A(j,17)>=0.04
        j=j+1;
    end
    while and(A(j,17)>0.03, A(j,17)<0.04);
        V0=V0+A(j,16); k=k+1; j=j+1;
    end
    Res(i,1) = V0./k; % MV0
    Vt=0; j=1; k=0; 
    while A(j,17)>=0.02  % here 0.02
        j=j+1;
    end
    l=j;
    while fl(j,2)==1
        Vt=Vt+A(j,16); k=k+1; j=j+1; 
    end
    Res(i,2) = Vt./k; % MVt
    
    tm=0; j=r-Term(i,1); k=0; 
    while A(j,16)>0.01 & k<Term(i,1)-1
        j=j+1;k=k+1;
    end
    if k==0
        k=1;
    end
    [Vm,tm] =min(A(end-Term(i,1):end-Term(i,1)+k-1, 16));
    if tm==1
        tm=2;
    end
    j=1; k=0; 
    while A(j,17)>=0.02
        j=j+1;
    end
    
    Res(i,3) = A(end-Term(i,1),16); % VT
    Res(i,4) = (Res(i,1)-Res(i,2))./Res(i,1); % speed decrease V0-Vt
    Res(i,5) = 1-(Res(i,2)-Res(i,3))./Res(i,2); % speed decrease Vt-VT 
    Res(i,6) = Vm;
    Res(i,7) = tm-1;
    Res(i,8) = (Res(i,3)-Res(i,6))./(Res(i,7)./1000); %a_l
    Res(i,9) = (mean(A(end-Term(i,1)-2.*fpf-3:end-Term(i,1)-2.*fpf,16))-A(end-Term(i,1),16))./(2.*fpf./1000);   % Aerodynamic deceleration
    Res(i,10) = Res(i,3)-Res(i,9).*(Res(i,7)./1000);
    Res(i,11) = A(end-Term(i,1)-fpf,17);                %distance from target at start of decceleration
    Res(i,12) = 180.*2*atan(0.00375./Res(i,11))./pi;     % Angular size at start of deceleration
    Res(i,13) = 0.00375./A(end-Term(i,1)-fpf,16);         %time to impact at start of deceleration
    Res(i,14) = A(lgt(i,1),17);
       if lgt(i,1)>19
    Res(i,15) = mean(A(lgt(i,1)-19:lgt(i,1),16));
    elseif lgt(i,1)<=19 & lgt(i,1)>9
    Res(i,15) = (A(lgt(i,1)-9:lgt(i,1),16));
    elseif lgt(i,1)<=9
        Res(i,15) = A(lgt(i,1),16);
    end
    
% figure(1)
% plot(A(:,1)-A(:,7),A(:,2)-A(:,8), 'k', 'LineWidth', 2.0)
% hold on; plot(0,0, 'or')
% axis([-0.15 0.15 -0.15 0.15])
% figure(2)
% plot(A(1:end-Term(i,1),17),A(1:end-Term(i,1),16),'.-','LineWidth', 2.0); 
% hold on
 t = (-(r-1-Term(i,1)).*0.001:0.001:Term(i,1).*0.001)';
%  figure(3) % velocity time
%  plot(t,A(:,16), '.-','LineWidth', 1.5);
%  hold on
    dt = [atan(0.00375./A(:,17)) atan(0.00375./A(:,17)) atan(0.00375./A(:,17))];
% figure(4) %subtence angle
% plot(t,dt(:,1), 'LineWidth', 1.5);
% hold on
dtt = velocity(dt,1000);
% figure(5) %angular speed
% plot(t,dtt(:,1), 'LineWidth', 1.5);
% hold on
cont = A(:,16).*sin(2.*dt(:,1))./(2.*A(:,17));
% figure(6) %angular velocity
% plot(1000.*A(:,17),cont, 'LineWidth', 1.5)
% hold on
% figure(7) % time to impact
% plot(A(:,17),A(:,17)./A(:,16),'LineWidth', 1.5)
%scatter(t,A(:,17)./A(:,16),'ok')
% hold on
% figure(8)
%plot(A(1:end-Term(i,1),17),A(1:end-Term(i,1),12),'LineWidth', 2.0)
% plot(A(:,17),A(:,12),'LineWidth', 2.0)
% hold on
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% [pks,locs]=findpeaks(A(:,16));
% plot(A(:,16))
% hold on
% plot(locs,pks, 'or')
% [pks,locs]=findpeaks(-A(:,16));
% plot(locs,-pks, 'ob')
% plot(r-Term(i,1),A(end-Term(i,1),16), '*k')
% plot(lgt(i,1),A(lgt(i,1),16),'*k')
% filenm = strcat('Fig',num2str(i));
% axis([1 r 0 0.6])
% legend(filenm)
% saveas(gcf,filenm,'emf')
% close
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% plot(A(:,17),A(:,16))
% hold on
% plot(A(r-Term(i,1),17),A(r-Term(i,1),16),'*k')
% plot(A(lgt(i,1),17),A(lgt(i,1),16),'*k')
% filenm = strcat('Fig',num2str(i));
% axis([0 A(1,17) 0 0.6])
% legend(filenm)
% saveas(gcf,filenm,'emf')
% close
Res(i,16) = (atan(0.0015./Res(i,14)))./cont(lgt(i,1),1);
Res(i,17) = cont(lgt(i,1),1);
Res(i,18) = mean(A(end-Term(i,1)-6.*fpf:end-Term(i,1)-4.*fpf,16));
Res(i,19) = mean(A(end-Term(i,1)-8.*fpf:end-Term(i,1)-6.*fpf,16));
Res(i,9) = (Res(i,18)-Res(i,3))./(4.*fpf./1000);   % Aerodynamic deceleration
Res(i,20) = A(end-Term(i,1)-3.*fpf,17);
Res(i,21) = A(1,17);
Res(i,22) = cont(lgt(i,1),1)./2.*dt(lgt(i,1),1);
Res(i,23) = dt(lgt(i,1),1)./cont(lgt(i,1),1);
Res(i,24) = 180.*cont(lgt(i,1),1)./pi;

    % for m=1:8
    %     Res(i,21+m) = A(end-Term(i,1)-m.*fpf,16);
    % end
% testing models
myfit = fittype('(x^(a+1))/b');
X=A(end-Term(i,1)-60:end-Term(i,1),17);
Y=A(end-Term(i,1)-60:end-Term(i,1),16);
coefs=fit(X,Y,myfit);
Ps(i,1)=coefs.a;
Es(i,1)=coefs.b;
% figure(9)
% plot(coefs,X,Y)
% hold on
end