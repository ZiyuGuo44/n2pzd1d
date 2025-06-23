% clear all
% 
% %load 1d_dz10_mld30_w1.mat
% %load 1d_dz10_mld30_w5.mat
% %load 1d_dz10_mld30_w10.mat
% %load 1d_dz10_mld30_w20.mat
% %load 1d_dz10_mld30_w30.mat
% 
% load 1d_dz5_mld30_w1.mat
% 
% close all

zdepths = deltaz * [1:ndepths];
[grid_days,grid_z] = meshgrid([1:ndays],zdepths);

cmap = jet;
myXlims = [0 (ndays/deltaday)];
myXtick = [0:(90/deltaday):(ndays/deltaday)];
myXtickLabel = [0:90:ndays];
    
myYlims = [0 (zmax/deltaz)];
myYtick = [0:(40/deltaz):(zmax/deltaz)];
myYtickLabel = [0:40:zmax];
    
for i =1:numExperiments
    mldSSP = plotData(i).mldSSP;  
    mld_layers = mldSSP/deltaz;  
    TodeSSP = plotData(i).TodeSSP; 
    kzSSP = plotData(i).kzSSP;
    parSSP = plotData(i).parSSP;

    NnodeSSP = plotData(i).NnodeSSP;
    NpodeSSP = plotData(i).NpodeSSP; 
    NcodeSSP = plotData(i).NcodeSSP; 
    P1nodeSSP = squeeze(plotData(i).PnodeSSP(1:ndepths,:));  % P1: [ndepths × ndays]
    P2nodeSSP = squeeze(plotData(i).PnodeSSP(ndepths+1:end,:));  % P2: [ndepths × ndays]
    ZnodeSSP = plotData(i).ZnodeSSP;
    DnodeSSP = plotData(i).DnodeSSP; 
    DcodeSSP = plotData(i).DcodeSSP; 
    DOMnodeSSP = plotData(i).DOMnodeSSP; 
    DOMcodeSSP = plotData(i).DOMcodeSSP; 

    pf_weighted1_tSSP = plotData(i).pf_weighted1_tSSP;
    pf_weighted2_tSSP = plotData(i).pf_weighted2_tSSP;
    
    C2P_phyto1_tSSP = plotData(i).C2P_phyto1_tSSP;
    C2P_phyto2_tSSP = plotData(i).C2P_phyto2_tSSP;
    C2P_detritus_tSSP = plotData(i).C2P_detritus_tSSP;
    C2P_DOM_tSSP = plotData(i).C2P_DOM_tSSP;
    C2P_zoo_tSSP = plotData(i).C2P_zoo_tSSP;
    
    GPPn1_tSSP = plotData(i).GPPn1_tSSP;
    GPPn2_tSSP = plotData(i).GPPn2_tSSP;

    GPPc1_tSSP = plotData(i).GPPc1_tSSP;
    GPPc2_tSSP = plotData(i).GPPc2_tSSP;
    GSPc1_adj_tSSP = plotData(i).GSPc1_adj_tSSP;    % homeostasis flux
    GSPc2_adj_tSSP = plotData(i).GSPc2_adj_tSSP;
    GSPc1_old_tSSP = plotData(i).GSPc1_old_tSSP;
    GSPc2_old_tSSP = plotData(i).GSPc2_old_tSSP; 
    GSPc1_new_tSSP = plotData(i).GSPc1_new_tSSP; 
    GSPc2_new_tSSP = plotData(i).GSPc2_new_tSSP; 

    POCfluxSSP = plotData(i).POCfluxSSP; 
    DOCfluxSSP = plotData(i).DOCfluxSSP; 

    interval=12;

 %========================================================================
 % physical parameters====================================================

    if i==1
        figure(1)

        subplot(1,3,1)
        contourf(TodeSSP,linspace(10,24,interval))
        set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
        set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
        axis square, axis ij
        grid on
        colormap(cmap)
        colorbar
        hold on
        plot (mld_layers,'w-','LineWidth',1);
        hold off
%        xlabel('Days','FontSize',14)
        xlabel('Days')
        title('','(a) Temperature (C)')
%        ylabel('Depth','FontSize', 14)
        ylabel('Depth')

        subplot(1,3,2)
        contourf(kzSSP,linspace(0,300,interval))
        set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
        set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
        axis square, axis ij
        colormap(cmap)
        grid on
        colorbar
        hold on
        plot (mld_layers,'w-','LineWidth',1);
        hold off
        title('','(b) KZ (m2/day)')
        xlabel('Days')

        subplot(1,3,3)
        contourf(parSSP,linspace(0,200,interval))
        set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
        set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
        axis square, axis ij
        grid on
        colormap(cmap)
        colorbar
        hold on
        plot (mld_layers,'w-','LineWidth',1);
        hold off
        title('','(c) PAR (W/m2)')
        xlabel('Days')
    end

    %return

 %========================================================================
 % stock=============================================================
    figure(2)

    min_val=0; 
    max_val=2.9;
    max_val1=3.4;
    interval=12;

    subplot(3,6,i)
    contourf(P1nodeSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ',num2str(i)],'P1 (mmolN*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(3,6,6+i)
    contourf(P2nodeSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    colormap(cmap)
    clim([min_val max_val]);
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','P2 (mmolN*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i == 6  
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.45, 0.02, 0.4]); 
    end

    subplot(3,6,12+i)
    contourf(P1nodeSSP+P2nodeSSP,linspace(min_val,max_val1,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val1]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','P1+P2 (mmolN*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i == 6  
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.1, 0.02, 0.25]); 
    end

 %========================================================================
 % preference=============================================================
    figure(3)

    min_val=0; 
    max_val=1;

    subplot(2,6,i)
    contourf(pf_weighted1_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ',num2str(i)],'P1 preference')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(2,6,6+i)
    contourf(pf_weighted2_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','P2 preference')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i == 6  
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.3, 0.02, 0.4]); 
    end

 %========================================================================
 % C:P=============================================================
    figure(4)

    min_val=80; 
    max_val=340;

    subplot(5,6,i)
    contourf(C2P_phyto1_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)], ' P1 C:P') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(5,6,6+i)
    contourf(C2P_phyto2_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','P2 C:P') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(5,6,12+i)
    contourf(C2P_detritus_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','D C:P') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(5,6,18+i)
    contourf(C2P_DOM_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','DOM C:P') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(5,6,24+i)
    contourf(C2P_zoo_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','Z C:P') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i == 6  
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.2, 0.02, 0.6]); 
    end

 %========================================================================
 % GSP_adj=============================================================
    figure(5)

    subplot(3,6,i)
    min_val=0; 
    max_val=0.65;

    contourf(ZnodeSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)], 'Z (mmolN*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.7, 0.02, 0.2])
    end

    subplot(3,6,6+i)
    min_val=0; 
    max_val=2.6;

    contourf(GSPc1_adj_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','Homeo1 adj (mmolC/m3/d)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end

    subplot(3,6,12+i)
    contourf(GSPc2_adj_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square
    axis ij
    clim([min_val max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','Homeo2 adj (mmolC/m3/d)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.15, 0.02, 0.4])
    end

 %========================================================================
 % GPP====================================================================
    figure(6)

    min_val=0; 
    max_val=0.6;

    subplot(4,6,i)
    contourf(GPPn1_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([min_val max_val]);
    %  colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off    
    title(['Experiment ', num2str(i)],'GPP1 (mmolN/m3/d)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end

    subplot(4,6,6+i)
    contourf(GPPn2_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([min_val max_val]);
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','GPP2 (mmolN/m3/d)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb1 = colorbar;
        set(cb1, 'Position', [0.92, 0.6, 0.02, 0.3])
    end

    min_val=0; 
    max_val=10;

    subplot(4,6,12+i)
    contourf(GPPc1_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([min_val max_val]);
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','GPP1 (mmolC/m3/d)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end

    subplot(4,6,18+i)
    contourf(GPPc2_tSSP,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([min_val max_val]);
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','GPP2 (mmolC/m3/d)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb3 = colorbar;
        set(cb3, 'Position', [0.92, 0.15, 0.02, 0.3])
    end

%===========================================================================
% nutrient and carbon ======================================================
    figure(7);

    zexpand = 70;
    n_layers = zexpand / deltaz;
    myYlims1 = [0 n_layers];
    myYtick1 = [0:(10/deltaz):n_layers];
    myYtickLabel1 = 0:deltaz:zexpand;

    subplot(4,6,i)
    contourf(NnodeSSP,linspace(0,9,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 9])
    %colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)],'NO3 (mmolN/m3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.76, 0.02, 0.15])
    end

    subplot(4,6,6+i)
    contourf(NnodeSSP(1:n_layers,:),linspace(0,0.8,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims1,'Ytick',myYtick1,'YtickLabel',myYtickLabel1)
    axis square, axis ij
    colormap(cmap)
    clim([0 0.8])
%    colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','NO3 (mmolN/m3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb1 = colorbar;
        set(cb1, 'Position', [0.92, 0.54, 0.02, 0.15])
    end

    subplot(4,6,12+i)
    contourf(NpodeSSP,linspace(0,0.5,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 0.5])
%    colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','PO4(mmolP/m3)') 
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.33, 0.02, 0.15])
    end

    subplot(4,6,18+i)
    contourf(NcodeSSP,linspace(1900,2050,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    clim([1900, 2050]);
    colormap(cmap)
%    colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','DIC (mmolC*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb3 = colorbar;
        set(cb3, 'Position', [0.92, 0.11, 0.02, 0.15])
    end

 %========================================================================
 % POC (stock and flux)===================================================
    figure(8)

    xmax = 17;

    subplot(2,6,i)
    contourf(DcodeSSP,linspace(0,xmax,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 xmax])
     % colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)],'POC (mmolC*m-3)')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.63, 0.02, 0.25])
        %set(cb, 'Position', [0.92, 0.4, 0.02, 0.25])
    end

    subplot(2,6,6+i)
    contourf(POCfluxSSP,linspace(0,xmax,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 xmax])
    %colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','POC flux(mmolC/m2/day)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.16, 0.02, 0.25])
        %set(cb2, 'Position', [0.92, 0.1, 0.02, 0.25])
    end

 %=========================================================================
 % DOC (stock and flux )===================================================
    figure(9);

%    xmax = 24; % for dz=10
    xmax = 22;  % for dz=5
    xmin = -4; % for dz=10, 5
%    xmax1 = 34; % for dz=10
    xmax1 = 35; % for dz=5
    xmax2 = 4; % for dz=10, 5

    subplot(3,6,i)
    contourf(DOMcodeSSP,linspace(0,xmax,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 xmax])
    %colorbar    
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)],'DOC (mmolC*m-3)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb = colorbar;
        set(cb, 'Position', [0.92, 0.73, 0.02, 0.2])
    end

    subplot(3,6,6+i)
    contourf(DOCfluxSSP,linspace(xmin,xmax1,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([0 xmax1])
    %colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','DOC flux(mmolC/m2/day)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.41, 0.02, 0.2])
    end

    subplot(3,6,12+i)
    contourf(DOCfluxSSP,linspace(xmin,xmax2,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    colormap(cmap)
    clim([xmin xmax2])
    %colorbar
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','DOC flux(mmolC/m2/day)')
    if i == 1  
        ylabel('Depth','FontSize', 14)
    end
    if i==6
        cb2 = colorbar;
        set(cb2, 'Position', [0.92, 0.12, 0.02, 0.2])
    end

 %=========================================================================
 % DOCflux profile ========================================================
    figure(10)

    day1 = 1;
    day2 = 100;
    day3 = 180;
    xmin = -4;
%    xmax = 34;  % for dz=10
    xmax = 36;

    subplot(2,3,1), hold on
        plot(POCfluxSSP(:,day1),zdepths)
    subplot(2,3,2), hold on
        plot(POCfluxSSP(:,day2),zdepths)
    subplot(2,3,3), hold on
        plot(POCfluxSSP(:,day3),zdepths)
    subplot(2,3,4), hold on
        plot(DOCfluxSSP(:,day1),zdepths)
    subplot(2,3,5), hold on
        plot(DOCfluxSSP(:,day2),zdepths)
    subplot(2,3,6), hold on
        plot(DOCfluxSSP(:,day3),zdepths)

    if i==6
        subplot(2,3,1)
            title('Day 1 POC flux')
            plot([xmin xmax], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(2,3,2)
            title('Day 100 POC flux')
            plot([xmin xmax], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(2,3,3)
            title('Day 180 POC flux')
            plot([xmin xmax], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(2,3,4)
            title('Day 1 DOC flux')
            plot([xmin xmax], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
       subplot(2,3,5)
            title('Day 100 DOC flux')
            plot([xmin xmax], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(2,3,6)
            title('Day 180 DOC flux')
            plot([xmin xmax], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m2/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
            legend('Exp 1','Exp 2,','Exp 3','Exp 4','Exp 5','Exp 6','MLD')
    end

%  GSP adjustment ratio=============================================================
    %  %========================================================================
    
    GSPc1_homeo2old = GSPc1_adj_tSSP./GSPc1_old_tSSP;
    GSPc2_homeo2old = GSPc2_adj_tSSP./GSPc2_old_tSSP;
    GSPct_homeo2old = (GSPc1_adj_tSSP+GSPc2_adj_tSSP)./(GSPc1_old_tSSP+GSPc2_old_tSSP);
    
    min_val = 0;
    max_val = 0.65;
    interval = 12;

    figure(11)

    subplot(3,6,i)
    contourf(GSPc1_homeo2old,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    clim([min_val, max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title(['Experiment ', num2str(i)],'GSPc1 homeo/total')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    subplot(3,6,6+i)
    contourf(GSPc2_homeo2old,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    clim([min_val, max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','GSPc2 homeo/total')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    
    subplot(3,6,12+i)
    contourf(GSPct_homeo2old,linspace(min_val,max_val,interval))
    set(gca,'Xlim',myXlims,'Xtick',myXtick,'XtickLabel',myXtickLabel)
    set(gca,'Ylim',myYlims,'Ytick',myYtick,'YtickLabel',myYtickLabel)
    axis square, axis ij
    clim([min_val, max_val]);
    colormap(cmap)
    hold on
    plot (mld_layers,'w-','LineWidth',1);
    hold off
    title('','GSPct homeo/total')
    if i == 1 
        ylabel('Depth','FontSize', 14)
    end
    
    if i == 6  
    cb = colorbar;
    set(cb, 'Position', [0.92, 0.3, 0.02, 0.4]); 
    end

 %=========================================================================
 % GSP profiles ========================================================
    
    GPPct = GPPc1_tSSP + GPPc2_tSSP;
    GSPct_adj = GSPc1_adj_tSSP + GSPc2_adj_tSSP;
    GSPct_old = GSPc1_old_tSSP + GSPc2_old_tSSP;
    GSPct_new = GSPc1_new_tSSP + GSPc2_new_tSSP;
    
    day1 = 1;
    day2 = 100;
    day3 = 180;
    xmin = 0;
    xmax = 10;
    xmax1 = 5.5;

    figure(12)

    subplot(4,3,1),  hold on, plot(GPPct(:,day1),zdepths)
    subplot(4,3,2),  hold on, plot(GPPct(:,day2),zdepths)
    subplot(4,3,3),  hold on, plot(GPPct(:,day3),zdepths)
    subplot(4,3,4),  hold on, plot(GSPct_old(:,day1),zdepths)
    subplot(4,3,5),  hold on, plot(GSPct_old(:,day2),zdepths)
    subplot(4,3,6),  hold on, plot(GSPct_old(:,day3),zdepths)
    subplot(4,3,7),  hold on, plot(GSPct_adj(:,day1),zdepths)
    subplot(4,3,8),  hold on, plot(GSPct_adj(:,day2),zdepths)
    subplot(4,3,9),  hold on, plot(GSPct_adj(:,day3),zdepths)
    subplot(4,3,10), hold on, plot(GSPct_new(:,day1),zdepths)
    subplot(4,3,11), hold on, plot(GSPct_new(:,day2),zdepths)
    subplot(4,3,12), hold on, plot(GSPct_new(:,day3),zdepths)

    if i==6
        subplot(4,3,1)
            title('Day 1 GPP')
            plot([xmin xmax], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,2)
            title('Day 100 GPP')
            plot([xmin xmax], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,3)
            title('Day 180 GPP')
            plot([xmin xmax], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,4)
            title('Day 1 GSP before homeostasis')
            plot([xmin xmax1], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,5)
            title('Day 100 GSP before homeostasis')
            plot([xmin xmax1], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,6)
            title('Day 180 GsP before homeostasis')
            plot([xmin xmax1], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,7)
            title('Day 1 Homeostasis')
            plot([xmin xmax1], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,8)
            title('Day 100 Homeostasis')
            plot([xmin xmax1], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,9)
            title('Day 180 Homeostasis')
            plot([xmin xmax1], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,10)
            title('Day 1 GSP after homeostasis')
            plot([xmin xmax1], [mldSSP(day1) mldSSP(day1)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,11)
            title('Day 100 GSP after homeostasis')
            plot([xmin xmax1], [mldSSP(day2) mldSSP(day2)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
        subplot(4,3,12)
            title('Day 180 GSP after homeostasis')
            plot([xmin xmax1], [mldSSP(day3) mldSSP(day3)],':')
            xlabel('mmolC/m3/day'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, axis square, box on, hold off
            legend('Exp 1','Exp 2,','Exp 3','Exp 4','Exp 5','Exp 6','MLD')
    end

 %=========================================================================
 % Stock partition profiles ===============================================
    
    TOTn  = NnodeSSP + P1nodeSSP + P2nodeSSP + ZnodeSSP + DnodeSSP + DOMnodeSSP;
    f_N   = NnodeSSP./TOTn;
    f_P1  = P1nodeSSP./TOTn;
    f_P2  = P2nodeSSP./TOTn;
    f_Z   = ZnodeSSP./TOTn;
    f_D   = DnodeSSP./TOTn;
    f_DOM = DOMnodeSSP./TOTn;

    day1 = 1;
    day2 = 180;
    xmin = 0;
    xmax = 1.1;
    xmax1 = 11;

    figure(13)

    subplot(3,7,1), hold on, plot(mean(TOTn,2),zdepths)
    subplot(3,7,2), hold on, plot(mean(f_N,2),zdepths)
    subplot(3,7,3), hold on, plot(mean(f_P1,2),zdepths)
    subplot(3,7,4), hold on, plot(mean(f_P2,2),zdepths)
    subplot(3,7,5), hold on, plot(mean(f_Z,2),zdepths)
    subplot(3,7,6), hold on, plot(mean(f_D,2),zdepths)
    subplot(3,7,7), hold on, plot(mean(f_DOM,2),zdepths)
    
    subplot(3,7,8), hold on, plot(TOTn(:,day1),zdepths)
    subplot(3,7,9), hold on, plot(f_N(:,day1),zdepths)
    subplot(3,7,10), hold on, plot(f_P1(:,day1),zdepths)
    subplot(3,7,11), hold on, plot(f_P2(:,day1),zdepths)
    subplot(3,7,12), hold on, plot(f_Z(:,day1),zdepths)
    subplot(3,7,13), hold on, plot(f_D(:,day1),zdepths)
    subplot(3,7,14), hold on, plot(f_DOM(:,day1),zdepths)
    
    subplot(3,7,15), hold on, plot(TOTn(:,day2),zdepths)
    subplot(3,7,16), hold on, plot(f_N(:,day2),zdepths)
    subplot(3,7,17), hold on, plot(f_P1(:,day2),zdepths)
    subplot(3,7,18), hold on, plot(f_P2(:,day2),zdepths)
    subplot(3,7,19), hold on, plot(f_Z(:,day2),zdepths)
    subplot(3,7,20), hold on, plot(f_D(:,day2),zdepths)
    subplot(3,7,21), hold on, plot(f_DOM(:,day2),zdepths)
        
    if i==6
        subplot(3,7,1)
            title('Total N','Annual mean')
            xlabel('mmolN/m3'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, box on, hold off
        subplot(3,7,2)
            title('NO3/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(3,7,3)
            title('P1/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(3,7,4)
            title('P2/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(3,7,5)
            title('Z/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(3,7,6)
            title('Detritus/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        subplot(3,7,7)
            title('DOM/Total N','Annual mean')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
            legend('Exp 1','Exp 2,','Exp 3','Exp 4','Exp 5','Exp 6')

        subplot(3,7,8)
            title('','Day 1')
            xlabel('mmolN/m3'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, box on, hold off
        for j=9:14
            subplot(3,7,j)
            title('','Day 1')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        end

        subplot(3,7,15)
            title('','Day 180')
            xlabel('mmolN/m3'), ylabel('Depth')
            axis([xmin xmax1 0 200]), axis ij, box on, hold off
        for j=16:21
            subplot(3,7,j)
            title('','Day 180')
            xlabel('Fraction'), ylabel('Depth')
            axis([xmin xmax 0 200]), axis ij, box on, hold off
        end
    end


end 

return