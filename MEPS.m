 % MEM Beamformer
function [w_meps] = MEPS(N,Rss,as,Theta_1,Theta,R_auto,Rs_auto)  
% tic
       u=[1,zeros(1,N-1)]';
        for nn=1:length(Theta)
            Theta_sv = exp(1j*pi * [0:N-1]'*sin(Theta(nn))); 
            nom=Theta_sv*Theta_sv';
            deno=(abs(Theta_sv'*inv(Rss)*u))^2;
            R_auto=R_auto+(nom/deno);
        end
        for nnn=1:length(Theta_1)
            A_Theta = exp(1j*pi * [0:N-1]'*sin(Theta_1(nnn)));
            nom1=A_Theta*A_Theta';
%             deno1=A_Theta'*inv(Rss)*(u*u')*inv(Rss')*A_Theta;
             deno1=(abs(A_Theta'*inv(Rss)*u))^2;
            Rs_auto=Rs_auto+(nom1/deno1);
        end
        ae_auto=sqrt(N)*Rs_auto*as;       % 13 Feb. 2020 sqrt(N)*
        w_meps=(1/real(ae_auto'*inv(R_auto)*ae_auto))*inv(R_auto) * ae_auto;
% toc