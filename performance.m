function [Power, SNR, BER] = performance(QPSK_g, N, QPSK_Re, QPSK_Im,ii,u_n1,sigma2_vec)


 %pre-allocation of matrix/vector
    QPSK_out = zeros(N/2,2);
    error = zeros(1,N/2);


    for i = 1: N/2

        d1 = (QPSK_g(i,1)-1)^2 +(QPSK_g(i,2)-1)^2; %1+j
        d2 = (QPSK_g(i,1)+1)^2 +(QPSK_g(i,2)-1)^2; %-1+j
        d3 = (QPSK_g(i,1)+1)^2 +(QPSK_g(i,2)+1)^2; %-1-j
        d4 = (QPSK_g(i,1)-1)^2 +(QPSK_g(i,2)+1)^2; %1-j

        [~,min_pos] = min([d1 d2 d3 d4]);

        switch min_pos
            case 1
                QPSK_out(i,1) = 1;
                QPSK_out(i,2) = 1;
            case 2
                QPSK_out(i,1) = -1;
                QPSK_out(i,2) = 1;
            case 3
                QPSK_out(i,1) = -1;
                QPSK_out(i,2) = -1;
            case 4
                QPSK_out(i,1) = 1;
                QPSK_out(i,2) = -1;


        end

        error(i) = (real(QPSK_Re(i))~=QPSK_out(i,1) || imag(QPSK_Im(i))~=QPSK_out(i,2))';

    end

    %BER of the transmitted signal
    BER = nnz(error)/length(error);
    %Power of the signal: (root-mean-squared)^2 of the signal
    signal = [real(u_n1),imag(u_n1)];
    Power = rms(signal)^2;
    SNR = Power/sigma2_vec(ii);

end