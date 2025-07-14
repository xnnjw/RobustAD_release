clearvars;
addpath '../SSR_algorithms/'
addpath '..'

%%
% % simulation parameters
M = 2^8; % = 128 (nr of atoms) 
K = 2^2; % = 4 sparsity 
L = 2^4; % = 16 (nr of snapshots) 
% dims = [2^3,2^4,2^5,2^6]; % Different choices of N
dims = [2^3]; % Different choices of N
SNR = 0:10;

nSNR = length(SNR);
meanrec = cell(1,length(dims));
meanTP = cell(1,length(dims));
meanTime = cell(1,length(dims));
GaussianS = true; 
% gaussianS = true if non-zero sources are Gaussian distributed
NRSIM = 2000;
for i  = 1:length(dims)
    N = dims(i);
    [meanrec{i},meanTP{i},meanTime{i}] = fnc_SNR_vs_PER(N,M,L,K,SNR,NRSIM,GaussianS);
end
% fname=sprintf('mat/SNRvsPER_N=%dto%d_K=%d_M=%d_L=%d_NRSIM=%d_GaussianS=%d_%s.mat',dims(1),dims(end),K,M,L,NRSIM,GaussianS,date);
% save(fname)
%%
plotData = true;
if plotData  
    nSNR = length(SNR);
    sigmas = zeros(nSNR,K);
    % Generate the source signal powers 2nd, 3rd, and 4th have SNR that is
    % lower than 1db, 2dB and 4dB than the first source
    for ii = 1:(K/4)
        sigmas(:,(ii-1)*4+1) = 10.^(SNR.'/10);
        sigmas(:,(ii-1)*4+2) = 10.^((SNR.'-1)/10);
        sigmas(:,(ii-1)*4+3) = 10.^((SNR.'-2)/10);
        sigmas(:,(ii-1)*4+4) = 10.^((SNR.'-4)/10);
    end    
    gm = geomean(sigmas,2);
    
    %% THIS IS THE PER PLOT 
    col2=[0.30100,0.74500,0.93300];
    col1=[0.46600,0.67400,0.18800];
    x = SNR;

    fignro=40;
    figure(fignro); clf
    subplot(1,4,1);
    hold on;
    plot(x,meanrec{1}(3,:),'r<-','DisplayName',sprintf('CL-OMP %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);    
    plot(x,meanrec{1}(4,:),'bv-','DisplayName',sprintf('CL-BCD %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{1}(5,:),'m^-','DisplayName',sprintf('CL-CCD %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{1}(2,:),'x-.','DisplayName',sprintf('SOMP %d',dims(1)),'LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanrec{1}(1,:),'ks-.','DisplayName',sprintf('SNIHT %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{1}(6,:),'c*-.','DisplayName',sprintf('CWOpt %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    legend('Location','northwest','FontSize',18);
    ylabel('Exact recovery')
    xlabel('SNR')
    grid on
    ylim([0,1])
       
    subplot(1,4,2)
    hold on;
    plot(x,meanrec{2}(3,:),'r<-','DisplayName',sprintf('CL-OMP %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{2}(4,:),'bv-','DisplayName',sprintf('CL-BCD %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{2}(5,:),'m^-','DisplayName',sprintf('CL-CCD %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{2}(2,:),'x--','DisplayName',sprintf('SOMP %d',dims(2)),'LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanrec{2}(1,:),'ks--','DisplayName',sprintf('SNIHT %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{2}(6,:),'c*-.','DisplayName',sprintf('CWOpt %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    %legend('Location','southeast','FontSize',18);
    xlabel('SNR')
    grid on
    ylim([0,1])
    
    subplot(1,4,3)
    hold on;
    plot(x,meanrec{3}(3,:),'r<-','DisplayName','CL-OMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{3}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{3}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{3}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanrec{3}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{3}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    xlabel('SNR')
    grid on
    ylim([0,1])

    subplot(1,4,4)
    hold on;
    plot(x,meanrec{4}(3,:),'r<-','DisplayName','CLOMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{4}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{4}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{4}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanrec{4}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanrec{4}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    %ylabel('Exact recovery')
    xlabel('SNR')
    grid on
    ylim([0,1])


    %% THIS IS THE TP PLOT 
    col2=[0.30100,0.74500,0.93300];
    col1=[0.46600,0.67400,0.18800];
    x = SNR;

    fignro=41;
    figure(fignro); clf
    subplot(1,4,1);
    hold on;
    plot(x,meanTP{1}(3,:),'r<-','DisplayName',sprintf('CL-OMP %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);    
    plot(x,meanTP{1}(4,:),'bv-','DisplayName',sprintf('CL-BCD %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{1}(5,:),'m^-','DisplayName',sprintf('CL-CCD %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{1}(2,:),'x-.','DisplayName',sprintf('SOMP %d',dims(1)),'LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTP{1}(1,:),'ks-.','DisplayName',sprintf('SNIHT %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{1}(6,:),'c*-.','DisplayName',sprintf('CWopt %d',dims(1)),'LineWidth',0.8,'MarkerSize',12);
    legend('Location','northwest','FontSize',18);
    ylabel('Exact recovery')
    xlabel('SNR')
    grid on
    ylim([0,1])
       
    subplot(1,4,2)
    hold on;
    plot(x,meanTP{2}(3,:),'r<-','DisplayName',sprintf('CL-OMP %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{2}(4,:),'bv-','DisplayName',sprintf('CL-BCD %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{2}(5,:),'m^-','DisplayName',sprintf('CL-CCD %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{2}(2,:),'x--','DisplayName',sprintf('SOMP %d',dims(2)),'LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTP{2}(1,:),'ks--','DisplayName',sprintf('SNIHT %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{2}(6,:),'c*-.','DisplayName',sprintf('CWopt %d',dims(2)),'LineWidth',0.8,'MarkerSize',12);
    %legend('Location','southeast','FontSize',18);
    xlabel('SNR')
    grid on
    ylim([0,1])
    
    subplot(1,4,3)
    hold on;
    plot(x,meanTP{3}(3,:),'r<-','DisplayName','CL-OMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{3}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{3}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{3}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTP{3}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{3}(6,:),'c*-.','DisplayName','CWopt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    % xlabel('SNR')
    grid on
    ylim([0,1])

    subplot(1,4,4)
    hold on;
    plot(x,meanTP{4}(3,:),'r<-','DisplayName','CLOMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{4}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{4}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{4}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTP{4}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTP{4}(6,:),'c*-.','DisplayName','CWopt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    %ylabel('Exact recovery')
    xlabel('SNR')
    grid on
    ylim([0,1])

    %% THIS IS THE CPU PLOT 
    col2=[0.30100,0.74500,0.93300];
    col1=[0.46600,0.67400,0.18800];
    x = SNR;

    fignro=42;
    figure(fignro); clf
    subplot(1,4,1);
    hold on;
    plot(x,meanTime{1}(3,:),'r<-','DisplayName','CL-OMP','LineWidth',0.8,'MarkerSize',12);    
    plot(x,meanTime{1}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{1}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{1}(2,:),'x-.','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTime{1}(1,:),'ks-.','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{1}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    legend('Location','northwest','FontSize',18);
    ylabel('Exact recovery')
    xlabel('SNR')
    title(sprintf('N = %d',dims(1)),'FontSize',18);
    grid on

    subplot(1,4,2)
    hold on;
    plot(x,meanTime{2}(3,:),'r<-','DisplayName','CL-OMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{2}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{1}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{2}(2,:),'x-.','DisplayName','OMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTime{2}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{2}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','southeast','FontSize',18);
    xlabel('SNR')
    title(sprintf('N = %d',dims(2)),'FontSize',18);
    grid on
    
    subplot(1,4,3)
    hold on;
    plot(x,meanTime{3}(3,:),'r<-','DisplayName','CL-OMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{3}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{3}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{3}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTime{3}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{3}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    xlabel('SNR')
    title(sprintf('N = %d',dims(3)),'FontSize',18);
    grid on

    subplot(1,4,4)
    hold on;
    plot(x,meanTime{4}(3,:),'r<-','DisplayName','CLOMP','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{4}(4,:),'bv-','DisplayName','CL-BCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{4}(5,:),'m^-','DisplayName','CL-CCD','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{4}(2,:),'x--','DisplayName','SOMP','LineWidth',0.8,'MarkerSize',12,'Color',col1);
    plot(x,meanTime{4}(1,:),'ks--','DisplayName','SNIHT','LineWidth',0.8,'MarkerSize',12);
    plot(x,meanTime{4}(6,:),'c*-.','DisplayName','CWOpt','LineWidth',0.8,'MarkerSize',12);
    %legend('Location','northwest','FontSize',18);
    %ylabel('Exact recovery')
    xlabel('SNR')
    title(sprintf('N = %d',dims(4)),'FontSize',18);
    grid on


end





