 clear;
clc;
rng(0);
n_sampl=30;                % Number of the snapshots
N=20;                      % Number of sensors
dl=0.5; % d/lambda
sqrt05 = sqrt(0.5);
% tet= [30 50];             % Location of interference users
tet= [-50  30];             % assumed Location of interference users
tet_s=10;                  % Presumed Location of desired user
tet_s_r=10;                 % Actual Location of desired user
m  = length(tet);          % Number of inteference users
p  = 1000;                 % Power of the inteference users

signal=[tet_s tet];    %Signals
NS=length(signal);  % Number of Signals

SNR_dB=-30:10:30;
for i=1:length(SNR_dB), SNR(i) = 10^(SNR_dB(i)/10); end
as= exp(1i*2.0*[0:N-1]'*pi*dl*sin(tet_s*pi/180.0));
as_r= exp(1i*2.0*[0:N-1]'*pi*dl*sin(tet_s_r*pi/180.0));
Am = exp(1i*2.0*[0:N-1]'*pi*dl*sin(tet*pi/180.0));           % inter. steering matr.

%-------------------------------------------------------------------------
% Quadratic constraint design 
for iter=1:100
    disp(iter);

wd=12;
wdd=4;
Theta=[-90*pi/180:0.9*pi/180:(tet_s-(wd/2))*pi/180    (tet_s+(wd/2))*pi/180:0.9*pi/180:90*pi/180 ];
Theta_1=(tet_s-(wd/2))*pi/180:0.9*pi/180:(tet_s+(wd/2))*pi/180;
%%%------------------------------------------------------------------------

   Ncount=0;     
  for ps=SNR
        ps;
        Ncount=Ncount+1;
        nk=n_sampl;
             
        Rss = zeros(N);
        Rx_k1=zeros(1);
        Rx_k2=zeros(1);
        Rs_auto=zeros(N);
        R_auto=zeros(N);
        zeta1=0;
        zeta11=0;
       
        for l=1:nk
            s = (randn(m,1) + 1i*randn(m,1)).*sqrt(p)'*sqrt05; % waveform
            ssig=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt05; 
            z = (randn(N,1) + 1i*randn(N,1))*sqrt05;           % noise sig^2=1
            xs = Am*s + ssig*as_r + z; % snapshot with signal (passive location case)
            zeta1=zeta1+(1/(n_sampl^2))*norm(xs)^4;
            zeta11=zeta11+(norm(xs))^4;
            Rss = Rss + xs*xs';
        end   
        
        for l=1:nk
            s_k = (randn(1,1) + 1i*randn(1,1)).*sqrt(p)'*sqrt05; % waveform
            ssig=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt05; 
            z = (randn(N,1) + 1i*randn(N,1))*sqrt05;           % noise sig^2=1
            xs_k1 = Am(:,1)'*Am(:,1)*s_k + Am(:,1)'*z; % snapshot with signal (passive location case)
            Rx_k1 = Rx_k1 + xs_k1*xs_k1';
        end   
        for l=1:nk
            s_k = (randn(1,1) + 1i*randn(1,1)).*sqrt(p)'*sqrt05; % waveform
            ssig=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt05; 
            z = (randn(N,1) + 1i*randn(N,1))*sqrt05;           % noise sig^2=1
            xs_k2 = Am(:,2)'*Am(:,2)*s_k + Am(:,2)'*z; % snapshot with signal (passive location case)
            Rx_k2 = Rx_k2 + xs_k2*xs_k2';
        end  
        
        Rss=Rss/nk;  % sample cov. matrix        
        R=p*(Am*Am')+eye(N,N);        % exact cov. matrix
%==================================================================
% MEPS
% tic
  [w_meps] = MEPS(N,Rss,as,Theta_1,Theta,R_auto,Rs_auto);
% toc
 Theta_22=(tet(1)-(wd/2))*pi/180:0.2*pi/180:(tet(1)+(wd/2))*pi/180; 
        thint_k=[];
        for l=1:nk
            Pr=[];
            for k=1:length(Theta_22) 
                aa=exp(1j*pi * [0:N-1]'*sin(Theta_22(k)));
                s = (randn(m,1) + 1i*randn(m,1)).*sqrt(p)'*sqrt05; % waveform
                ssig=(randn(1,1) + 1i*randn(1,1))*sqrt(ps)*sqrt05; 
                z = (randn(N,1) + 1i*randn(N,1))*sqrt05;           % noise sig^2=1
                xs = Am*s + ssig*as_r + z; % snapshot with signal (passive location case)
                pr=abs(xs'*aa);
                Pr=[Pr pr];
            end
            [Prmax, I]=max(Pr);
            thint=Theta_22(I)*180/pi;
            thint_k=[thint_k thint];
            delta_tet=2*nk*(thint_k-tet(1))/((2*l)-nk);
        end  
        
        Delta=max(ceil(delta_tet));
        max(thint_k);
%==============================================================================
%%% Proposed Method 
 Sum_C=0; 
 sum=0;
 [Es D_s]=eig(Rss);
        for ii=NS+1:N
            sum=sum+D_s(ii,ii);
        end
 sig_n=(1/(N-NS))*sum;  %average Noise power
 gama1=(1/((N-NS)*2*pi))*sum;
 sig_int1=(Rx_k1- sig_n*(Am(:,1)'*Am(:,1)))/(abs(Am(:,1)'*Am(:,1))^2); % Pwr of int 1
 sig_int2=(Rx_k2- sig_n*(Am(:,2)'*Am(:,2)))/(abs(Am(:,2)'*Am(:,2))^2); % Pwr of int 2
 Intsig=[sig_int1  sig_int2];
 gama2=max(Intsig); 
 
%  Theta_2c=(tet(1)-(2))*pi/180:0.9*pi/180:(tet(1)+(2))*pi/180;
 Theta_2c=min(thint_k)*pi/180:0.9*pi/180:max(thint_k)*pi/180;
 Theta_3c=(tet(2)-(2))*pi/180:0.9*pi/180:(tet(2)+(2))*pi/180;
 Theta_ic=[Theta_2c Theta_3c];

 % maximum of interference powers
        for nc=1:length(Theta_ic)
            Thetai_C = exp(1j*pi * [0:N-1]'*sin(Theta_ic(nc))); 
            sum_c=Thetai_C*Thetai_C';
            Sum_C=Sum_C+sum_c;  % Interference covariance matrix
            C=2*pi*gama1*eye(N)+(gama2-gama1)*Sum_C;
        end
%  Ra=C; 
 w_c=(1/real(as'*inv(C)*as))*inv(C) * as;

%============================================================================        
        SINRopt(Ncount,iter)=ps*as_r'*inv(R)*as_r;                                      % Optimal
        SINR_c(Ncount,iter)=ps*(abs(w_c'*as_r)*abs(w_c'*as_r))/abs(w_c'*R*w_c); % proposed
        SINRmeps(Ncount,iter)=ps*(abs(w_meps'*as_r)*abs(w_meps'*as_r))/abs(w_meps'*R*w_meps); % New Steering vector
      
  end
end

SINRopt=mean(SINRopt,2);
SINR_c=mean(SINR_c,2);
SINRmeps=mean(SINRmeps,2);


snr=-30:10:30;
sinr_o=SINRopt.';
sinr_c=SINR_c.';
sinr_ME=SINRmeps.';


OPTMAL=10*log10(abs(sinr_o))
DRS=10*log10(abs(sinr_c))
MEPS=10*log10(abs(sinr_ME))

lw=1.5;
msz=8;
figure(1)

plot(snr,OPTMAL,'k',snr,DRS,'-b>',snr,MEPS,'-o','LineWidth',lw,'MarkerSize',msz);
axis([-30 30  -30 50])
xlabel('SNR (dB)','FontSize',14)
ylabel('SINR (dB)','FontSize',14)
legend('Optimal','IPN-PSEUR','IPN-MEPS')
grid on

%================================================PSD Plot
% fc = 1e6; %carrier freqiency
% lambda = 3*10^8/fc; %wavelength
% d = 0.5*lambda;%distance between neigbouring antennas
% 
% teta_v = -pi/2:pi/720:pi/2; %range of scanned directions in radians
% teta_v_deg = teta_v/pi*180; %same as above but in degrees
% Nteta = length(teta_v); %number of scanned directions
% 
% iiRxx=inv(Rtilde_in);
% iRxx=inv(R_auto);
% iRhat=inv(Rss);
% for k=1:Nteta  %scanning loop
%     teta = teta_v(k); %current scannig direction
%     delta_fi = 2*pi*d/lambda*sin(teta); %relateive phase delay 
%     for m=1:N
%         aT(m) = exp(1i*((m-1)*delta_fi)); %steering vector
%     end
%    
%      a=aT.'; %transpose of steering vector
%      
%       wMVDR = (iiRxx*a)/(a'*iiRxx*a); %weight vector for Caponbeamforming
%       PMVDR(k) = abs(1/(a'*iRhat*a));
%      
%       wAuto = (iRxx*a)/(a'*iRxx*a);
%       Pmeps(k) = abs(1/(a'*iRhat*(u*u')*iRhat'*a));
% end
% 
% PMVDR_dB = 10*log10(PMVDR);
% Pmeps_dB=  10*log10(Pmeps);
% figure (2)
% lw=1.5;
% plot(teta_v_deg, PMVDR_dB,'--r',teta_v_deg, Pmeps_dB,'b','LineWidth',lw)
% axis([-90 90 -30 60])
% % xticks([-90 -50  -20   3  15  50 90])
% % xticklabels({'-90','-50',' -20','3','15','50','90'})
% legend('Capon','MEM');
% % title('spatial power spectrum of Capon beamformer')
% set(gca,'FontSize',13);
% xlabel('\theta, Degree','FontSize',15);ylabel('Power Spectrum','rotation',90,'FontSize',15);
% 
% % axes('position',[.624 .40 .28 .28])
% % box on % put box around new pair of axes
% % indexOfInterest = (teta_v_deg < 18) & (teta_v_deg > 0);
% % plot(teta_v_deg(indexOfInterest),PMVDR_dB(indexOfInterest),'--r',teta_v_deg(indexOfInterest),Pauto_dB(indexOfInterest),'b','LineWidth',lw) % plot on new axes
% % set(gca,'FontSize',8);
% % axis tight
% grid on

