clear all;
num_itr = 100000;
ct = 100;
R = 0.5;
eps = 2^R-1;
T=0.5;
M=32;
PL = eps; PH = (1+PL)*eps;
snrdb = [0 : 5: 20];

for isnr =  1: length(snrdb)
    aoi_tdma_ana(isnr)=0;
    aoi_tdma_mp_ana(isnr)=0;
    aoi_noma_ana(isnr)=0;
    aoi_noma_ana_mp(isnr)=0;
    for m = 1 : M/2
        m=4;mp = m+M/2;
        snr=10^(snrdb(isnr)/10);snrS = snr;

        aoi_tdma_ana(isnr) =aoi_tdma_ana(isnr)+ m*T + M*T*(2*exp(eps/snr)-1)/2;
        aoi_tdma_mp_ana(isnr) = aoi_tdma_mp_ana(isnr)+ mp*T + M*T*(2*exp(eps/snr)-1)/2;


        %analytical results -  user m
        ptemp = exp(-eps/snrS)*(1-exp(-(snr*eps/snrS+1)*eps/snr))/(snr*eps/snrS+1);
        p0 = (1-exp(-eps/snr)-ptemp)*(1-exp(-eps/snrS)/(1+snr*eps/snrS))+...
            ptemp*(1-exp(-eps/snrS));
        pm = exp(-eps/snrS);
        pmp = (1-exp(-eps/snr)-ptemp)*(exp(-eps/snrS)/(1+snr*eps/snrS))+...
            ptemp*(exp(-eps/snrS));

        deltam0 = 1/(pm+pmp)^2*(1-p0)^2*((pm+pmp)*m*T*pm/(1-p0)^2   + pmp/2*m*T*pm/(1-p0)...
                                       + (pm+pmp)*mp*T*pmp/(1-p0)^2 - pm/2*mp*T*pmp/(1-p0));

        aoi_noma_ana(isnr) = aoi_noma_ana(isnr)+deltam0 +(M^2*T^2*(pm+pmp)^2*(1+p0)+M^2/2*T^2*pm*pmp*(1-p0)^2)...
            /M/T/(pm+pmp)^2/(1-p0)/2;

        %analytical results -  user m'
        p0 = (1-exp(-eps/snrS)/(1+eps*snr/snrS))*(1-exp(-eps/snr));
        pm = exp(-eps/snrS)/(1+eps*snr/snrS);
        pmp = (1-exp(-eps/snrS)/(1+eps*snr/snrS))*(exp(-eps/snr));

        deltam0 = 1/(pm+pmp)^2*(1-p0)^2*((pm+pmp)*m*T*pm/(1-p0)^2   + pmp/2*m*T*pm/(1-p0)...
                                       + (pm+pmp)*mp*T*pmp/(1-p0)^2 - pm/2*mp*T*pmp/(1-p0));

        aoi_noma_ana_mp(isnr) = aoi_noma_ana_mp(isnr)+ deltam0 +(M^2*T^2*(pm+pmp)^2*(1+p0)+M^2/2*T^2*pm*pmp*(1-p0)^2)...
            /M/T/(pm+pmp)^2/(1-p0)/2;
    end
    
end
%plot(snrdb,aoi_tdma,snrdb,aoi_tdma_ana,  snrdb,aoi_noma_m,snrdb,aoi_noma_ana)
plot(snrdb,(aoi_tdma_ana+aoi_tdma_mp_ana)/M, snrdb,(aoi_noma_ana+aoi_noma_ana_mp)/M)
 
