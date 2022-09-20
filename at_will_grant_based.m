clear all;
num_itr = 2000;
ct = 100;
R = 1;
eps = 2^R-1;
T=0.5;
M=8;
PL = eps; PH = (1+PL)*eps;
snrdb = [0 : 5: 20 ];

for isnr =  1: length(snrdb)
    
    snr=10^(snrdb(isnr)/10);snrS = snr;
    %TDMA - generate num_itr frames
    h = complex(sqrt(0.5)*randn(num_itr,1),sqrt(0.5)*randn(num_itr,1));   %m  
    h = abs(h).^2;
    h_sts = max(0,sign(h-eps/snr)); h_sts(1)=1; h_sts(end)=1; %the use of special caes at the two ends have no impact if num_itr is large 
    start_vec = strfind(h_sts',[1 0])'; % begining of one update
    end_vec = strfind(h_sts',[0 1])'; % ending of one update
    succ_vec = strfind(h_sts',[1 1])'; % these are continous success 
    xj =  end_vec - start_vec(1:length(end_vec)) + 1;
    xj = [xj; ones(length(succ_vec),1)];
    aoi_tdma(isnr) = sum(T*xj*M*T + (xj*M*T).^2/2)/sum(xj*M*T);

    aoi_tdma_ana(isnr) = T + M*T*(2*exp(eps/snr)-1)/2;
    aoi_tdma_app(isnr) = T + M*T/2;
     
    
    %NOMA- generate num_itr frames - user m
    h = complex(sqrt(0.5)*randn(num_itr,4),sqrt(0.5)*randn(num_itr,4));   %two users' channel gains in two time slots in each frame
    h = abs(h).^2;
    number_qj = 0;
    start_q = 0; % the calculation starts from the beginning
    Qj = []; yj = [0];
    ind_position = 1;% 1 means (j-1)-th update at m, 0 means (j-1)-th update at m'
    sum0=0; sum1=0;sum2=0;
    for i = 1 : num_itr
        %decide user m is a weak or strong user at the first time slot
        if log2(1+snr*h(i,1))>R%successful update in the first time slot 
            yj(end) = yj(end)+M*T/2 + ind_position * M*T/2; % if start m, MT; if start m', MT/2
            ind_position = 1;
            yj = [yj 0]; %start a new update 
            Qj = [Qj T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            number_qj = number_qj +1;
            sum1 = sum1+1;
        elseif  log2(1+snrS*h(i,3)/(1+snr*h(i,4)))>R %successful update in the first time slot 
            yj(end) = yj(end)+ind_position * M*T/2  +  M*T; % if start m, 3MT/2; if start m', MT
            ind_position = 0;
            yj = [yj 0]; %start a new update 
            Qj = [Qj T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj 
            number_qj = number_qj +1;
            sum2 = sum2+1;
        else
            yj(end) = yj(end)+M*T;
            sum0 = sum0+1;
        end
    end
    aoi_noma_m(isnr) = sum(Qj(2:end-1))/sum(yj(2:end-1));

    p0 = (1-exp(-eps/snr))*(1-exp(-eps/snrS)/(1+snr*eps/snrS)) ;
    pm = exp(-eps/snrS);
    pmp = (1-exp(-eps/snr))*(exp(-eps/snrS)/(1+snr*eps/snrS)) ;

    aoi_noma_ana(isnr) = T +(M^2*T^2*(pm+pmp)^2*(1+p0)+M^2/2*T^2*pm*pmp*(1-p0)^2)...
        /M/T/(pm+pmp)^2/(1-p0)/2;
    
    aoi_noma_app(isnr) = T + M*T/2;

end
%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana,snrdb,aoi_noma_old,snrdb, aoi_noma_old_ana ,snrdb,aoi_noma_m,snrdb,aoi_noma_ana)
%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana,snrdb,aoi_noma_m,snrdb,aoi_noma_ana)
plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana, snrdb,aoi_tdma_app,snrdb,aoi_noma_m,snrdb,aoi_noma_ana,snrdb,aoi_noma_app)
%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana,snrdb,aoi_noma_m,snrdb,aoi_noma_ana)
 