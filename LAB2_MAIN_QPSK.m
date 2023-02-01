% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% LAB. EXERCISE 2: MODELING AND SIMULATION OF A SATELLITE LINK (2022/23)

% Students: David Pérez López & Daniel Talavera Jiménez

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Generation of N random bits
clc; clear;close all

N = 8000;
bits = rand(1,N)>0.5;

%% Build constellation

QPSK = bits*2-1;
QPSK_Re = QPSK(1:2:end-1);
QPSK_Im = QPSK(2:2:end)*1j;

%% Small scale fading

alpha = 0.2;
BW = 24e6; %[Hz]
T = (1+alpha)/BW; %[s]
f = 2e9; %[Hz]
tau = [1e-7,2e-7];

QPSK = QPSK_Re + QPSK_Im;


% for ii = 1: length(tau)
%     for nn = 0:length(QPSK)-1
% 
%         t = nn*T;
% 
%         if t>0
% 
%             u_t = -t/T +1;
%             p(nn+1,ii) = u_t - 0.3*(-(t-tau(ii))/T +1);
% 
%         elseif t<0
% 
%             u_t = t/T -1;
%             p(nn+1,ii) = u_t - 0.3*(-(t-tau(ii))/T +1);
% 
%         else
%             u_t = 1;
%             p(nn+1,ii) = u_t - 0.3*(-(t-tau(ii))/T +1);
%         end
% 
%     end
% end
% 
% A_p1 = conv(QPSK,p(:,1));
% A_p2 = conv(QPSK,p(:,2));

p_func = @(n) 1-n*T/T - 0.3*(1-(n*T-tau(1))/T); % IT DOES NOT PROVIDE THE DESIRED COEFFICIENTS


numerator1 = [1 0 -0.3];
numerator2 = [1 0 0 0 -0.3];

A_p1 = filter(numerator1,1,QPSK); %filter([numerator],[denominator],QPSK)
A_p2 = filter(numerator2,1,QPSK); %filter([numerator],[denominator],QPSK)

figure()
hold on
plot(real(A_p1),imag(A_p1),'o',LineWidth=1.6)
xlabel('\Re')
ylabel('\Im')
plot(real(A_p2),imag(A_p2),'*',LineWidth=1)
plot(real(QPSK),imag(QPSK),'ko','MarkerSize',15)
title('A[n]*p[n] for QPSK constellation')
legend('\tau_1=1e-7','\tau_2 =2e-7','Original constellation')
grid on


%% Add noise

% sigma2_vec = [10 5 1 0.5 0.1 0.01 0.001]; %array of standard deviations of the noise (noise according to a Gaussina distribution)
sigma2_vec = [linspace(1e-5,10), linspace(10.5,100,10)];

for ii = 1: length(sigma2_vec)

    sigma2 = sigma2_vec(ii); %amplitude of the noise

    z = sqrt(sigma2)*(randn(1,N)+1j*randn(1,N))/sqrt(2); % this needs to be added to the vector of symbols

    check_sigma2 = sum(abs(z).^2)/N; %this should be equal to sigma2

    % q = A*p+z;
    q_Re = real(A_p1)+real(z(1:2:end-1));
    q_Im = imag(A_p1)+imag(z(2:2:end));
    q_Re2 = real(A_p2)+real(z(1:2:end-1));
    q_Im2 = imag(A_p2)+imag(z(2:2:end));


if ii == 1
        figure()
        plot(QPSK_Re,imag(QPSK_Im),'r*',LineWidth=1)
        xlabel('$\Re$','Interpreter','latex');
        ylabel('$\Im$','Interpreter','latex');
        title(['QPSK constellation with \sigma^2 = ',num2str(0)])
        xlim([-1.2 1.2])
        ylim([-1.2 1.2])
        grid on

end
    
%     figure()
%     plot(q_Re,q_Im,'o',LineWidth=2)
%     xlabel('$\Re$','Interpreter','latex');
%     ylabel('$\Im$','Interpreter','latex');
%     title(['\sigma^2 = ',num2str(sigma2_vec(ii)),'(\tau_1)'])
%     axis square
% 
%     figure()
%     plot(q_Re2,q_Im2,'o',LineWidth=2)
%     xlabel('$\Re$','Interpreter','latex');
%     ylabel('$\Im$','Interpreter','latex');
%     title(['\sigma^2 = ',num2str(sigma2_vec(ii)),'(\tau_2)'])
%     axis square
%     
   %% Equalizer w[n]



    NN = 2;
   LL = 3;
   LLL = 5;
   NNN = 4;

   p0_2 = numerator2(1);  p0_1 = numerator1(1);
   p1_2 = numerator2(2);  p1_1 = numerator1(2);
   p2_2 = numerator2(3);  p2_1 = numerator1(3);
   p3_2 = numerator2(4);
   p4_2 = numerator2(5);

   P_1 = [p0_1 0 0;... 
                p1_1 p0_1 0; ...
                p2_1 p1_1 p0_1;...
                0 p2_1 p1_1; ...
                0 0 p2_1];

%     P_2 = [p0_2 0 0 0 0;... 
%                 p1_2 p0_2 0 0 0; ...
%                 p2_2 p1_2 p0_2 0 0;...
%                 0 p2_2 p1_2 p0_2 0; ...
%                 0 0 p2_2 p1_2 p0_2;...
%                 0 0 0 p2_2 p1_2;...
%                 0 0 0 0 p2_2];
 P_2 = [p0_2 0 0 0 0;... 
                p1_2 p0_2 0 0 0; ...
                p2_2 p1_2 p0_2 0 0;...
                p3_2 p2_2 p1_2 p0_2 0; ...
                p4_2 p3_2 p2_2 p1_2 p0_2;...
                0    p4_2 p3_2 p2_2 p1_2;...
                0    0    p4_2 p3_2 p2_2;...
                0    0    0    p4_2 p3_2;...
                0    0    0    0    p4_2 ];

    c1 = [1; zeros(NN+LL-1,1)];
    c2 = [1; zeros(NNN+LLL-1,1)];

    w_ZF1 = pinv(P_1)*c1;
    w_ZF2 = pinv(P_2)*c2;




    q1 = q_Re+q_Im*1j;
    u_n1 = filter(w_ZF1,1,q1); %filter([numerator],[denominator],QPSK)

    q2 = q_Re2+q_Im2*1j;
    u_n2 = filter(w_ZF2,1,q2); %filter([numerator],[denominator],QPSK)


%     figure()
%     plot(real(u_n1),imag(u_n1),'o',LineWidth=2)
%     xlabel('$\Re$','Interpreter','latex');
%     ylabel('$\Im$','Interpreter','latex');
%     title(['u[n] with \tau_1 = ',num2str(tau(1)),' and \sigma^2 =',num2str(sigma2_vec(ii))])
%     axis square
% 
%     figure()
%     plot(real(u_n2),imag(u_n2),'o',LineWidth=2)
%     xlabel('$\Re$','Interpreter','latex');
%     ylabel('$\Im$','Interpreter','latex');
%     title(['u[n] with \tau_2 = ',num2str(tau(2)),' and \sigma^2 =',num2str(sigma2_vec(ii))])
%     axis square




    %% Convert signal with noise back to the symbols with the original constellation
    % For this, geometric distance to each respective point of the constellation must be
    % calculated

    %Merge QPSK into a 'geometric' vector w/o imaginary parts


    QPSK_g = [real(u_n1)',imag(u_n1)'];
    QPSK2_g = [real(u_n2)',imag(u_n2)'];
    QPSK3_g = [q_Re',q_Im'];
    QPSK4_g = [q_Re2',q_Im2'];


    [Power1_i, SNR1_i, BER1_i] = performance(QPSK_g, N, QPSK_Re, QPSK_Im,ii,u_n1,sigma2_vec);
    [Power2_i, SNR2_i, BER2_i] = performance(QPSK2_g, N, QPSK_Re, QPSK_Im,ii,u_n2,sigma2_vec);
    [Power3_i, SNR3_i, BER3_i] = performance(QPSK3_g, N, QPSK_Re, QPSK_Im,ii,q1,sigma2_vec);
    [Power4_i, SNR4_i, BER4_i] = performance(QPSK4_g, N, QPSK_Re, QPSK_Im,ii,q2,sigma2_vec);
    Power(ii,:) = [Power1_i Power2_i Power3_i Power4_i];
    SNR(ii,:) = [SNR1_i SNR2_i SNR3_i SNR4_i]; 
    BER(ii,:) = [BER1_i BER2_i BER3_i BER4_i];

end

%% Plot of BER vs SNR
figure()
subplot(2,1,1)
semilogy(10*log10(SNR(:,1)),BER(:,1),LineWidth=1.5)
hold on
semilogy(10*log10(SNR(:,3)),BER(:,3),LineWidth=1.5)
xlabel('SNR [dB]','Interpreter','latex')
ylabel('BER (log scale)','Interpreter','latex')
grid on
legend('w equalizer','w/o equalizer')
title('QPSK with \tau_1 ')


% figure()
subplot(2,1,2)
semilogy(10*log10(SNR(:,2)),BER(:,2),LineWidth=1.5)
hold on
semilogy(10*log10(SNR(:,4)),BER(:,4),LineWidth=1.5)
xlabel('SNR [dB]','Interpreter','latex')
ylabel('BER (log scale)','Interpreter','latex')
grid on
legend('w equalizer','w/o equalizer')
title('QPSK with \tau_2')

