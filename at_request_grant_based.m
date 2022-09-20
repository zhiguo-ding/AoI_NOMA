clear all;
num_itr = 5000;
ct = 100;
R = 1;
eps = 2^R-1;
T=0.5;
M=8;
PL = eps; PH = (1+PL)*eps;
snrdb = [0 : 5: 20];

for isnr =  1: length(snrdb)
    
    m=1;mp = m+M/2;
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
    aoi_tdma(isnr) = sum(m*T*xj*M*T + (xj*M*T).^2/2)/sum(xj*M*T);
    aoi_tdma_mp(isnr) = sum(mp*T*xj*M*T + (xj*M*T).^2/2)/sum(xj*M*T);

    aoi_tdma_ana(isnr) = m*T + M*T*(2*exp(eps/snr)-1)/2;
    aoi_tdma_mp_ana(isnr) = mp*T + M*T*(2*exp(eps/snr)-1)/2;
    
    aoi_tdma_app(isnr) = m*T+M*T/2;
    aoi_tdma_mp_app(isnr) = mp*T+M*T/2; 
     
        
        
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
            yj = [yj 0]; %start a new update 
            if ind_position==1 % the previous cae is m
                Qj = [Qj m*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            else
                Qj = [Qj mp*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            end
            number_qj = number_qj +1;
            sum1 = sum1+1;
            ind_position = 1;
        elseif (log2(1+snrS*h(i,2)/(1+snr*h(i,1)))<R & log2(1+snrS*h(i,3)/(1+snr*h(i,4)))>R) | (log2(1+snrS*h(i,2)/(1+snr*h(i,1)))>R & log2(1+snrS*h(i,3))>R)%successful update in the first time slot 
            yj(end) = yj(end)+ind_position * M*T/2  +  M*T; % if start m, 3MT/2; if start m', MT
            yj = [yj 0]; %start a new update 
            if ind_position==1 % the previous cae is m
                Qj = [Qj m*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj 
            else
                Qj = [Qj mp*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            end
            number_qj = number_qj +1;
            sum2 = sum2+1;
            ind_position = 0;
        else
            yj(end) = yj(end)+M*T;
            sum0 = sum0+1;
        end
    end
    aoi_noma_m(isnr) = sum(Qj(2:end-1))/sum(yj(2:end-1));

    %NOMA- generate num_itr frames - user m'
    number_qj = 0;
    start_q = 0; % the calculation starts from the beginning
    Qj = []; yj = [0];
    ind_position = 1;% 1 means (j-1)-th update at m, 0 means (j-1)-th update at m'
    sum0=0; sum1=0;sum2=0;
    for i = 1 : num_itr
        %decide user m is a weak or strong user at the m-th time slot
        if log2(1+snr*h(i,2)/(1+snr*h(i,1)))>R%successful update in the first time slot 
            yj(end) = yj(end)+M*T/2 + ind_position * M*T/2; % if start m, MT; if start m', MT/2
            yj = [yj 0]; %start a new update 
            if ind_position==1 % the previous cae is m
                Qj = [Qj m*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            else
                Qj = [Qj mp*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            end
            number_qj = number_qj +1;
            sum1 = sum1+1;
            ind_position = 1;
        elseif (log2(1+snrS*h(i,2)/(1+snr*h(i,1)))<R & log2(1+snrS*h(i,4))>R)  %successful update in the m' time slot 
            yj(end) = yj(end)+ind_position * M*T/2  +  M*T; % if start m, 3MT/2; if start m', MT
            yj = [yj 0]; %start a new update 
            if ind_position==1 % the previous cae is m
                Qj = [Qj m*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj 
            else
                Qj = [Qj mp*T*yj(end-1)+yj(end-1)^2/2]; %using the updated yj to find Qj
            end
            number_qj = number_qj +1;
            sum2 = sum2+1;
            ind_position = 0;
        else
            yj(end) = yj(end)+M*T;
            sum0 = sum0+1;
        end
    end
    aoi_noma_mp(isnr) = sum(Qj(2:end-1))/sum(yj(2:end-1));
    
    %analytical results -  user m
    ptemp = exp(-eps/snrS)*(1-exp(-(snr*eps/snrS+1)*eps/snr))/(snr*eps/snrS+1);
    p0 = (1-exp(-eps/snr)-ptemp)*(1-exp(-eps/snrS)/(1+snr*eps/snrS))+...
        ptemp*(1-exp(-eps/snrS));
    pm = exp(-eps/snrS);
    pmp = (1-exp(-eps/snr)-ptemp)*(exp(-eps/snrS)/(1+snr*eps/snrS))+...
        ptemp*(exp(-eps/snrS));
    
    deltam0 = 1/(pm+pmp)^2*(1-p0)^2*((pm+pmp)*m*T*pm/(1-p0)^2   + pmp/2*m*T*pm/(1-p0)...
                                   + (pm+pmp)*mp*T*pmp/(1-p0)^2 - pm/2*mp*T*pmp/(1-p0));

    aoi_noma_ana(isnr) = deltam0 +M*T*(2*(pm+pmp)^2*(1+p0)+pm*pmp*(1-p0)^2)...
        /(pm+pmp)^2/(1-p0)/4;
    aoi_noma_app(isnr) = m*T+M*T/2;
    
    %analytical results -  user m'
    p0 = (1-exp(-eps/snrS)/(1+eps*snr/snrS))*(1-exp(-eps/snr));
    pm = exp(-eps/snrS)/(1+eps*snr/snrS);
    pmp = (1-exp(-eps/snrS)/(1+eps*snr/snrS))*(exp(-eps/snr));
    
    deltam0 = 1/(pm+pmp)^2*(1-p0)^2*((pm+pmp)*m*T*pm/(1-p0)^2   + pmp/2*m*T*pm/(1-p0)...
                                   + (pm+pmp)*mp*T*pmp/(1-p0)^2 - pm/2*mp*T*pmp/(1-p0));

    aoi_noma_ana_mp(isnr) = deltam0 +M*T*(2*(pm+pmp)^2*(1+p0)+pm*pmp*(1-p0)^2)...
        /(pm+pmp)^2/(1-p0)/4;
    aoi_noma_app_mp(isnr) = (m*T+mp*T*eps)/(1+eps)+M*T/2;
                                      
end
%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana, snrdb, aoi_noma_old ,snrdb,aoi_noma_old_ana,snrdb,aoi_noma_m,snrdb,aoi_noma_ana)
plot(snrdb,aoi_tdma_mp,snrdb,aoi_tdma_mp_ana,snrdb,aoi_tdma_mp_app, ...
    snrdb,aoi_noma_mp,snrdb,aoi_noma_ana_mp,snrdb,aoi_noma_app_mp)

%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana,snrdb,aoi_tdma_app, '-o','MarkerSize',10);
%plot(snrdb,aoi_noma_m,snrdb,aoi_noma_ana,snrdb,aoi_noma_app, '-x','MarkerSize',10);

%m
%plot(snrdb,ones(1,length(snrdb))*mp*T+M*T/2)
%mp NOMA
%plot(snrdb,ones(1,length(snrdb))*(m*T+mp*T*eps)/(1+eps)+M*T/2)