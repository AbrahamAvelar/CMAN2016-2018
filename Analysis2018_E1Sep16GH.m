% del plato 1:8 y 9:16 están repetidos en competencia los 8 arreglos de la
% CMAN. Los platos 17 y 18 CREO que son el arreglo 1 del CMAN pero en cultivos
% individuales para analizarse por OD y Los 19 y 20 son el arreglo 2. 

%%
load Data

%%
CleanData=BgDataAll;
for pl=1:length(BgDataAll)
    for t=1:length(BgDataAll(pl).t)
        if length(CleanData(pl).t)>t %porque una vez que quita uno de CleanData puede regarla al final
            val=nanmean(CleanData(pl).OD(t,:));
            if val<0.001;
                CleanData(pl).OD(t,:) = [];
                CleanData(pl).t(t)    = [];
                if pl<3
                    CleanData(pl).YFP(t,:)= [];
                    CleanData(pl).RFP(t,:)= [];
                end
            end
        end
    end
end

%% %poner las horas correctamente
%horasInoculos='C:\Users\JAbraham\Dropbox\PhD\Experimentos (ExMonth15)\E1Sep16 Longevidad en KO de CMAN2.0 via de tor\HorasInoculosEMOct16.xlsx';
horasInoculos='C:\Users\JAbraham\Documents\PRE-DROPBOX\Experimentos (ExMonth15)\E1Sep16 Longevidad en KO de CMAN2.0 via de tor\HorasInoculosEMOct16.xlsx';

[NUM TXT]= xlsread(horasInoculos);
odTh=.3;%cu'anto tiene que bajar la od en promedio para que sea considerado otro dia
for pl=[1:20] %ESTA PARTE ACOMODA LAS HORAS DE INOCULACION
    
    NuevosDias=EncuentraDias(CleanData(pl), odTh);
    ind=find(NUM(:,1)==pl)';
    
    y=cell2mat(TXT(ind,2));
    dias=datenum([y(:,end-1:end) y(:,3) y(:,4:6) y(:,1:2)  ]);
    diahora=dias+NUM(ind,3)-.001;%hay que hacer que sea un poquito mas chica para que no haya dos mediciones con la misma hora
    
    CleanData(pl).t(NuevosDias)=diahora;
    
end

%% %%%% EN VEZ DE PONERLE LA HORA AL 1ER PLATO, PONDRÉ OTRO QUE SEA OD, CFP Y RFP = 0.1 
horasInoculos='C:\Users\JAbraham\Dropbox\PhD\Experimentos (ExMonth15)\E1Sep16 Longevidad en KO de CMAN2.0 via de tor\HorasInoculosEMOct16.xlsx';
[NUM TXT]= xlsread(horasInoculos);
odTh=.3;%cu'anto tiene que bajar la od en promedio para que sea considerado otro dia
for pl=[1:20] %ESTA PARTE ACOMODA LAS HORAS DE INOCULACION
    NuevosDias=EncuentraDias(CleanData(pl), odTh);
	ind=find(NUM(:,1)==pl)';
	y=cell2mat(TXT(ind,2));
    dias=datenum([y(:,end-1:end) y(:,3) y(:,4:6) y(:,1:2) ]) ;
	diahora=dias+NUM(ind,3);    
    for i=1:length(NuevosDias)
        x(1:96)=.05;
        if i==1
            CleanData(pl).OD=[x; CleanData(pl).OD];
            CleanData(pl).RFP=[x; CleanData(pl).RFP];
            CleanData(pl).YFP=[x; CleanData(pl).YFP];
            CleanData(pl).t=[diahora(i) CleanData(pl).t];
        else
            
            CleanData(pl).OD = [CleanData(pl).OD(1:NuevosDias(i)+i-2,:); x; CleanData(pl).OD(NuevosDias(i)+i-1:end,:)];
            CleanData(pl).RFP = [CleanData(pl).RFP(1:NuevosDias(i)+i-2,:); x; CleanData(pl).RFP(NuevosDias(i)+i-1:end,:)];
            CleanData(pl).YFP = [CleanData(pl).YFP(1:NuevosDias(i)+i-2,:); x; CleanData(pl).YFP(NuevosDias(i)+i-1:end,:)];
            CleanData(pl).t=[CleanData(pl).t(1:NuevosDias(i)+i-2) diahora(i) CleanData(pl).t(NuevosDias(i)+i-1:end)];
        end
    end
end
%% MUTS esta es el arreglo general de las mutantes en el ensayo completo
for hide = 1% Genera Muts con sus posiciones
Muts.Tpk1=[17 74 12 53 64 94];
Muts.Tpk2=[25 82 28 61 15 72];
Muts.Msn2=[33 3 36 69 23 80];
Muts.Msn4=[49 11 44 85 39];
Muts.Rim15=[57 19 60 93 47];
Muts.Ras2=[65 35 68 22 55];
Muts.Gis1=[10 43 76 30 71];
Muts.Tor1=[18 51 92 46 79 4];
Muts.Sch9=[26 67 5 54 87];
Muts.yak1=[42 75 21 62 32];
Muts.rtg1=[50 83 29 78 40];
Muts.bmh1=[58 37 14 86 48];
Muts.Refs=[2 6 9 13 20 24 27 31 34 38 41 45 52 56 59 63 66 70 73 77 84 88 91 95];
Muts.YFPOnly=[1 16 89 90];
Muts.RFPOnly=[7 8 81 96];

end
% STRAINS
strains={'S288c(Mosaic/Lab)','L-1528(Europe/Ferment)','YJM981(Europe/Clinical)','YIIc17-E5(Mosaic/Ferment)','YPS606(American/Wild)','YJM975(Europe/Clinical)','Y12(Asia/Sake)','Y55(African/Lab)','S288c(Mosaic/Lab)','L-1528(Europe/Ferment)','YJM981(Europe/Clinical)','YIIc17-E5(Mosaic/Ferment)','YPS606(American/Wild)','YJM975(Europe/Clinical)','Y12(Asia/Sake)','Y55(African/Lab)','CREO-S288c(Mosaic/Lab)','CREO-S288c(Mosaic/Lab)','CREO-L-1528(Europe/Ferment)','CREO-L-1528(Europe/Ferment)'};
% ODths
for hide=1;
ODTh(1,:)=[.2 .75];
ODTh(2,:)=[.2 .74];
ODTh(3,:)=[.2 .58];
ODTh(4,:)=[.2 .65];
ODTh(5,:)=[.2 .6];%.64
ODTh(6,:)=[.2 .5];
ODTh(7,:)=[.2 .6];
ODTh(8,:)=[.2 .6];
ODTh(9,:)=[.2 .75];
ODTh(10,:)=[.2 .74];
ODTh(11,:)=[.2 .58];
ODTh(12,:)=[.2 .65];
ODTh(13,:)=[.2 .6];%.64
ODTh(14,:)=[.2 .5];
ODTh(15,:)=[.2 .6];
ODTh(16,:)=[.2 .6];
ODTh(17,:)=[.2 .7];
ODTh(18,:)=[.2 .7];
ODTh(19,:)=[.13 .5];
ODTh(20,:)=[.13 .5];
end

%% DECAY FROM SERGIO's MASTER
pls=[1:20]; %number of plates
od = .28; % od at interpolation ORIGINAL .28
odTh = -0.3; %ORIGINAL -0.3
n0 = .18; % minimum OD for GR ORIGINAL .18
nt = .55; % Max OD for GR ORIGINAL .48
% 3. Bkg medium substraction
value = .07;
for i=1:20
    value = min([min(CleanData(i).OD) value]);
end

bgdataClean = bkgsubstractionOD(CleanData,value);% Generate new matrix bgdataAll.OD
[timeOD] = getimeAA(bgdataClean,pls,od,odTh);  %Gets matrix with time at od stated above for each well and plate
[gwrate] = getgr(bgdataClean,pls,n0,nt,odTh); %Gets growth rate matrix
[survival] = getsurv(timeOD,gwrate,pls); %Gets survival percentage matrix
survival = getxvec(bgdataClean,survival,pls,od,odTh); %Gets matrix into survival structure for time at which interpolation was made in days and plots percentage vs time
% [decayRate] = getExponentialRate(survival,1:20); %esta necesita CurveFit que no lo puedo usar bien en mi compu

%%  SURVIVAL CURVES, MEDIA Y ERROR
%close all
con=0;
names=fieldnames(Muts);
for i=1:length(names)
mut=Muts.(cell2mat(names(i)));
refs=Muts.Refs;
titulo=names(i);
figure(i+100)
clf
hold on    
con=0;
for pl=1:20%16 o hasta 20
    %for well=swr1
        con=con+1;
        subplot(5,4,con)
        %plot(survival(pl).s(:,well+1),'o-r')
        stderror=nanstd(survival(pl).s(:,mut)')./sqrt(length(survival(pl).s(:,mut)'));
        meansurv=nanmean(survival(pl).s(:,mut)');
        NuevosDias=EncuentraDias(CleanData(pl), .3);
        x=CleanData(pl).t(NuevosDias)-CleanData(pl).t(1);
        a(1)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.8 .4 .10]);
        hold on
        stderror=nanstd(survival(pl).s(:,refs)')./sqrt(length(survival(pl).s(:,refs)'));
        meansurv=nanmean(survival(pl).s(:,refs)');
        a(2)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.5 .5 .5]);
        ylim([20 110])
        %xlim([0 15])
        title(strcat(titulo, ' PL: ', num2str(pl), strains(pl) ))
        grid on
    %end
end %8x12 curvas de decaimiento

legend(a,'mut','ref')
end
%%  SURVIVAL CURVES, MEDIA Y ERROR solo de los platos con cultivo individual
%close all
con=0;
names=fieldnames(Muts);
names=fieldnames(MutsDos);
for i=1:length(names)
% mut=Muts.(cell2mat(names(i)));
% refs=Muts.Refs;
% titulo=names(i);
figure(i+100)
clf
hold on    
con=0;
for pl=17:18%20%16 o hasta 20
mut=MutsDos(1).(cell2mat(names(i)));
refs=MutsDos(1).Refs;
titulo=names(i);

    %for well=swr1
        con=con+1;
        subplot(2,2,con)
        %plot(survival(pl).s(:,well+1),'o-r')
        stderror=nanstd(survival(pl).s(:,mut)')./sqrt(length(survival(pl).s(:,mut)'));
        meansurv=nanmean(survival(pl).s(:,mut)');
        NuevosDias=EncuentraDias(CleanData(pl), .3);
        x=CleanData(pl).t(NuevosDias)-CleanData(pl).t(1);
        a(1)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.8 .4 .10]);
        hold on
        stderror=nanstd(survival(pl).s(:,refs)')./sqrt(length(survival(pl).s(:,refs)'));
        meansurv=nanmean(survival(pl).s(:,refs)');
        a(2)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.5 .5 .5]);
        plot(x, survival(pl).s(:,mut),'color',[.8 .4 .10] )
        
        ylim([20 110])
        %xlim([0 15])
        title(strcat(titulo, ' PL: ', num2str(pl), strains(pl) ))
        grid on
    %end
end %8x12 curvas de decaimiento

legend(a,'mut','ref')
end
%%  SURVIVAL CURVES, MEDIA Y ERROR solo de los platos con cultivo individual
%close all
con=0;
Muts.RefsOrilla=[6 24 41 56 73 91 ];
Muts.RefsEsquina=[2   9   88  95];
Muts.RefsCentro=[ 13 20  27 31 34 38 45 52 59 63 66 70 77 84];

Muts.YFPOnlyNoEsquina=[ 16 90];
Muts.RFPOnlyNoEsquina=[7 81 ];names=fieldnames(Muts);
for i=18:length(names)
mut=Muts.(cell2mat(names(i)));
refs=Muts.Refs; %esto va normalmente
refs=[Muts.RefsEsquina Muts.RefsOrilla]; %esto na-mas pa comparar con las esquinas de los Onlys
titulo=names(i);
figure(i+100)
clf
hold on    
con=0;
for pl=17:20%16 o hasta 20
        con=con+1;
        subplot(2,2,con)
        stderror=nanstd(survival(pl).s(:,mut)')./sqrt(length(survival(pl).s(:,mut)'));
        meansurv=nanmean(survival(pl).s(:,mut)');
        NuevosDias=EncuentraDias(CleanData(pl), .3);
        x=CleanData(pl).t(NuevosDias)-CleanData(pl).t(1);
        a(1)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.8 .4 .10]);
        hold on
        stderror=nanstd(survival(pl).s(:,refs)')./sqrt(length(survival(pl).s(:,refs)'));
        meansurv=nanmean(survival(pl).s(:,refs)');
        a(2)=PlotConError(x,meansurv, meansurv+stderror,meansurv-stderror,[.5 .5 .5]);
        ylim([20 110])
        title(strcat(titulo, ' PL: ', num2str(pl), strains(pl) ))
        grid on
end %8x12 curvas de decaimiento

legend(a,'mut','ref')
end
%% ODBgDataAll y TBgDataAll CompareSCoeffs

odTh=.3;
ExOD=.35; %OD de extrapolacion
toremove=[Muts.RFPOnly Muts.YFPOnly];
previousS4=0;
refs=Muts.Refs;
tfijo = 9;
figuras=0;
horasExponencialEsp=[5.5 10; 5.4 9.5; 7.25 11.5; 6.2 9.75; 6.25 9.5; 5.5 9.5; 5.5 9.5; 5.3 9.5; 6 10.5; 5.5 9.25; 7.5 10.75; 6.5 10.5; 6.75 10.25; 6.25 10; 5.75 9.5; 5.75 10];
tfijos=[8.5 8 10 8.5 8.5 8.5 8 8 9 7.8 9.5 9 8 8.8 8.5 8.5];
%horasExponencial=[6 11];
%FluorMut='RFP';
%FluorRef='YFP';
pozos=[1:10];
for platos=1:16
    minOD=ODTh(platos,1);
    maxOD=ODTh(platos,2);
    CleanData(platos).CFP=CleanData(platos).YFP;
    tfijo=tfijos(platos);
    horasExponencial=horasExponencialEsp(platos,:);
    
	BgDataAll = calculaTiempos(CleanData, platos, odTh);
    
    temp=RoCtFijo(BgDataAll,tfijo, horasExponencial, odTh, platos, pozos, figuras);
    temp= CalcRelatSurv(temp, platos,odTh,pozos,platos);
    TBgDataAll(platos) =temp(platos);
    
    temp = RoCODFija(BgDataAll, ExOD, minOD, maxOD, platos,odTh);
    temp = CalcRelatSurv(temp, platos,odTh,pozos,platos+100);
    ODBgDataAll(platos)=temp(platos);
    %[TBgDataAll(platos), ODBgDataAll(platos) NSBgDataAll(platos)] = CompareSCoeffs(CleanData, odTh,tfijo, ExOD, minOD, maxOD, pozos, platos, horasExponencial, figuras, toremove,refs, previousS4)
    
end
close all
%%

platos=1;
toremove=[1 8 9];
refs=0;
col='g';
pozos=1:20;
A = CalcNSModel(ODBgDataAll,platos,pozos,toremove,refs,col);

%%
%%%%GENERAR EL HashDos MutsDos
[NUM TXT RAW]=xlsread('Copia_PosicionesCepas_Curado.xlsx');
for i=1:length(TXT)
    MutsDos(NUM(i,1)).(cell2mat((TXT(i)))) = NUM(i, 2:length(find(~isnan(NUM(i, 1:end)))));
	MutsDos(NUM(i,1)+8).(cell2mat((TXT(i)))) = NUM(i, 2:length(find(~isnan(NUM(i, 1:end)))));
end

%%
%save 2Mar18
load 2Mar18

%% FIGURA PARA IR DE .RoC A .s CON UNA SOLA MUTANTE EN EL FOR INTERNO
%platos=[3 2 6 4 1 5];%solo para la figura de tpk1/ras2/tpk2/yak1
DATA=TBgDataAll;
names=fieldnames(MutsDos);
platos=1;%[3 2 6 4 1]; % platos=[3 6 2 1 4 5];
con=0;
con2=0;
figure(233)
clf
cepas={'Ras2'};
for j=1:length(cepas)
    cepa=cepas(j);
for pl=platos%(6)
    con=con+1;
    subplot(length(cepas), length(platos),con)
    hold off
    useit=[2 3 4 6 7 12 16];
    mref=nanmedian(DATA(pl).s(Muts.Refs(useit),2));
    refs=DATA(pl).s(Muts.Refs( useit),2);
    plot([0 15], [zeros(length(refs),1) 15*(refs-mref)],'k');
    hold on
    mutantes = MutsDos(pl).(str2mat(cepa));%MutsDos(pl).Ras2;% Muts.Ras2
    for w = mutantes
        m=nanmean([DATA(pl).s(w,2)])%-mref;
        plot([0 1 15], [0 m m*15],'-r')
        plot(DATA(pl).Tdays(1:size(DATA(pl).RoC(:,w))), DATA(pl).RoC(:,w)-DATA(pl).s(w,1) ,'or')
    end %hace plot de las lineas de cada pozo
    for hide=1 %genera toplots 
        temp=nanmean([DATA(pl).s(mutantes,2)'; DATA(pl+8).s(mutantes,2)'])%-mref ;
        toplot2(con,1:length(temp))=temp;
        toplotStd(con)=nanstd(temp)/sqrt(length(temp));
        %temp=nanmean([CleanData(pl).s(useit,2)'; CleanData(pl+8).s(useit,2)'])-mref ;
        temp=refs-mref;
        toplot2R(con,1:length(temp))=temp;
        toplotStdR(con)=nanstd(temp)/sqrt(length(temp));
    end
    axis square
    %title( strcat(strains(pl), '-', cepa ))
    title( strcat(strains(pl)))
    ylabel('ln(RFP/YFP)')
    xlabel('Days')
    ylim([-1.5 1.5])
end
end
legend(cepas)

%% FIGURA PARA IR DE .RoC A .s CON UNA SOLA MUTANTE EN EL FOR INTERNO
%platos=[3 2 6 4 1 5];%solo para la figura de tpk1/ras2/tpk2/yak1
DATA=TBgDataAll;
names=fieldnames(MutsDos);
platos=[3 2 6 4 1]; % platos=[3 6 2 1 4 5];
con=0;
con2=0;
figure(233)
clf
cepas={'Tpk1','Tpk2','Ras2'};
cepas={'Tor1'};
for j=1:length(cepas)
    cepa=cepas(j);
for pl=platos%(6)
    con=con+1;
    subplot(length(cepas), length(platos),con)
    hold off
    useit=[2 3 4 6 7 12 16];
    mref=nanmedian(DATA(pl).s(Muts.Refs(useit),2));
    refs=DATA(pl).s(Muts.Refs( useit),2);
    plot([0 15], [zeros(length(refs),1) 15*(refs-mref)],'k');
    hold on
    mutantes = MutsDos(pl).(str2mat(cepa));%MutsDos(pl).Ras2;% Muts.Ras2
    for w = mutantes
        %m=nanmean([DATA(pl).s(w,2) DATA(pl+8).s(w,2)])-mref;
        m=nanmean([DATA(pl).s(w,2)])-mref;
        plot([0 1 15], [0 m m*15],'-r')
    end %hace plot de las lineas de cada pozo
    for hide=1 %genera toplots 
        temp=nanmean([DATA(pl).s(mutantes,2)'; DATA(pl+8).s(mutantes,2)'])-mref ;
        toplot2(con,1:length(temp))=temp;
        toplotStd(con)=nanstd(temp)/sqrt(length(temp));
        %temp=nanmean([CleanData(pl).s(useit,2)'; CleanData(pl+8).s(useit,2)'])-mref ;
        temp=refs-mref;
        toplot2R(con,1:length(temp))=temp;
        toplotStdR(con)=nanstd(temp)/sqrt(length(temp));
    end
    axis square
    %title( strcat(strains(pl), '-', cepa ))
    title( strcat(strains(pl)))
    ylabel('ln(RFP/YFP)')
    xlabel('Days')
    ylim([-1.5 1.5])
end
end
legend(cepas)

%% cacula una matriz con medias y medianas de 12platos*15kos
campos=fieldnames(MutsDos)
DATA=TBgDataAll;
platos=[1:6 9:14];
for i = 1:length(campos)
    for pl=platos
        muts=MutsDos(pl).(str2mat(campos(i)));
        sMean(i,pl) = nanmean(DATA(pl).s(muts,2));
        sMedian(i,pl)=nanmedian(DATA(pl).s(muts,2));
        srefs=nanmean(ODBgDataAll(pl).s(MutsDos(pl).Refs,2));
        NormsMean(i,pl) = nanmean(DATA(pl).s(muts,2))-srefs;
        NormsMedian(i,pl) = nanmedian(DATA(pl).s(muts,2))-srefs;
    end
end%%
for hide=1
% esto tal vez est[a mal porque imagesc voltea las matrices
% figure(64)
% clf
% imagesc(NormsMedian)
% title('Medianas Normalizadas')
% xlabel('plato')
% set(gca,'ytick',1:15,'yticklabel',campos)
% figure(67)
% clf
% imagesc(NormsMean)
% title('Media Normalizadas')
% xlabel('plato')
% set(gca,'ytick',1:15,'yticklabel',campos)
end
%% barras por GEN
platos=[1 9 2 10 3 11 4 12 6 14]; %acomodados por par de replica
figure(33)
clf
con1=0;
for i = [1 2 6 10 8 9 5 11 3 4 7 12];%1:12%13%length(campos)
    con1=con1+1;
    subplot(3,4,con1)
    con=0;
    for pl=platos
        con=con+1;
        muts=MutsDos(pl).(str2mat(campos(i)));
        srefs=nanmean(DATA(pl).s(MutsDos(pl).Refs,2));
        NormsMean(i,pl) = nanmean(DATA(pl).s(muts,2))-srefs;
        NormsMedian(i,pl) = nanmedian(DATA(pl).s(muts,2))-srefs;
        bar(con,NormsMean(i,pl),'g')
        hold on
        plot(con,(DATA(pl).s(muts,2))-srefs,'.r')
        if ~mod(con1,4)
         text(con,-.03,strains(pl),'Rotation',270)
        end
    end
    title(campos(i))
    ylim([-.08614 .08614])
    xlim([0 length(platos)+1])
end%%

%% barras por CEPA
platos=[1 9 2 10 3 11 4 12 5 13 6 14]; %acomodados por par de replica
%platos=[1 9 2 10 3 11 6 14 7 15 8 16]; %acomodados por par de replica
mutantes=[1 2 6 10 8 9 5 11 3 4 7 12];
figure(43)
clf
con1=0;
for pl=platos %1:12%13%length(campos)
    con1=con1+1;
    subplot(3,4,con1)
    %figure(con1)
    con=0;
    for i = mutantes
        con=con+1;
        muts=MutsDos(pl).(str2mat(campos(i)));
        srefs=nanmean(DATA(pl).s(MutsDos(pl).Refs,2));
        NormsMean(i,pl) = nanmean(DATA(pl).s(muts,2))-srefs;
        NormsMedian(i,pl) = nanmedian(DATA(pl).s(muts,2))-srefs;
        bar(con,NormsMean(i,pl),'g')
        hold on
        plot(con,(DATA(pl).s(muts,2))-srefs,'.r')
        if con1>8%~mod(con1,4)
         text(con,-.08614,campos(i),'Rotation',270)
        end
    end
    title(strains(pl))
    ylim([-.08614 .08614])
    xlim([0 length(mutantes)+1])
    set(gca, 'xtick',[])
    
end%%


%% 12 NotBoxPlot. uno por plato
%platos=[1 9 2 10 3 11 4 12 5 13 6 14]; %acomodados por par de replica
mutantes=[1 2 6 10 8 9 5 11 3 4 7 12];
figure(43)
clf
con1=0;
for pl=platos %1:12%13%length(campos)
    con1=con1+1;
    subplot(3,4,con1)
    %figure(con1)
    con=0;
    for i = mutantes
        con=con+1;
        muts=MutsDos(pl).(str2mat(campos(i)));
        srefs=nanmean(DATA(pl).s(MutsDos(pl).Refs,2));
        NormsMean(i,pl) = nanmean(DATA(pl).s(muts,2))-srefs;
        NormsMedian(i,pl) = nanmedian(DATA(pl).s(muts,2))-srefs;
        %bar(con,NormsMean(i,pl),'g')
        hold on
        %plot(con,(DATA(pl).s(muts,2))-srefs,'.r')
        NotBoxPlotAA((DATA(pl).s(muts,2))-srefs,ones(1,length(muts))*con)
        hold on
        %if ~mod(con1,4)
         text(con,-.03,campos(i),'Rotation',270)
        %end
    end
    title(strains(pl))
    ylim([-.08614 .08614])
    xlim([0 length(mutantes)+1])
end%%
%% En un solo BoxPlot Todas juntas por cepa
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(34)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6 13];
platos=[3 2 6 4 1 5]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.45:.075:.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=colormap;
DATA=TBgDataAll;
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2

        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        refere=DATA(pl).s(MutsDos(pl).Refs,2);
        srefs=nanmean(refere);
        toplot(con,1:length(ind))=(DATA(pl).s(ind,2))-srefs;
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(i)-.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        hold on
        text(ones(1,length(ind))*con+diferencias(i)-.013, toplot(con,:), {ind} )
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.01
            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end
        
        
        
        ind = MutsDos(pl+8).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);
        refere=DATA(pl+8).s(MutsDos(pl+8).Refs,2);
        srefs=nanmean(refere);
        toplot2(con,1:length(ind+8))=(DATA(pl+8).s(ind,2))-srefs;
        toplot2(toplot==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        hold on        
        text(ones(1,length(ind))*con+diferencias(i)+.013, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.01
            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('slope,s')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.152 .13090])
%ylim([-.12 .085])
xlim([0.5 5.5])% Solo 5 cepas 
xlim([0.5 6.5])
legend(names(orden),'Location', 'Best')

%% En un solo BoxPlot Todas juntas por gen
figure(489)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1]; % platos=[3 6 2 1 4 5];
diferencias=-.375:.18:.45; %diferencias=-.375:.15:.45;
mapita=colormap(HSV);
DATA=TBgDataAll;

for j=1:6;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*10,:), 'MarkerFaceColor', mapita(j*10,:) );
    hold on
end %Es para que la leyenda quede bien
 con=0;
for i = 1:length(orden)
    con=con+1;
    con2=0;
    for pl=platos;
        clear toplot toplot2
        con2=con2+1;
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        
        refere=DATA(pl).s(MutsDos(pl).Refs,2);
        srefs=nanmean(refere);
        toplot(con,1:length(ind))=(DATA(pl).s(ind,2))-srefs;
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(con2)-.03,0.001,'patch',20,[1,1,1],mapita(con2*10,:) ,[.3 .3 .3] ,3)
        hold on
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.05
            plot( con+diferencias(con2)-.03, max(toplot(con,:))+.01, '*r')
        end
        
        refere=DATA(pl+8).s(MutsDos(pl+8).Refs,2);
        srefs=nanmean(refere);
        toplot2(con,1:length(ind))=(DATA(pl+8).s(ind,2))-srefs;
        toplot2(toplot2==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(con2)+.03,0.001,'patch',20,[1,1,1],mapita(con2*10,:) ,[.3 .3 .3] ,3)
        hold on
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.05
            plot( con+diferencias(con2)+.03, max(toplot2(con,:))+.01, '*r')
        end
    end
    hold on
    plot( [-1 13], [0 0] )
    set(gca, 'Xtick', 1:12, 'xticklabel', [])
	%h=text( 1:12, zeros(12,1)-.151,names(orden) );
    %set( h, 'rotation', 270 )
    set(gca, 'XtickLabel',names(orden) );
    ylabel('slope,s')
end
for i=1:13
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.75 .75 .75] )
end
ylim([-.172 .1175])
xlim([0.5 12.5])
legend(strains(platos),'Location', 'Best')
title('por gen')
%%
% solo los de la parte de ras
xlim([8.5 12.5])
ylim([-.16 .085])

%% En un solo BoxPlot Tres/Dos cepas interesantes
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(332)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6 13];%[5 9 8 2 1 6 13]%
platos=[1 2 5]%[3 2 5]Casos extremos de posible sesgo%[2 1]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.45:.075:.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=colormap;
DATA=TBgDataAll;
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        refere=DATA(pl).s(MutsDos(pl).Refs,2);
        srefs=nanmean(refere);
        toplot(con,1:length(ind))=(DATA(pl).s(ind,2))-srefs;
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(i)-.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.9,.9,.9],7)
        hold on
        text(ones(1,length(ind))*con+diferencias(i)-.031, toplot(con,:), {ind} )

        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.01
            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end

        ind = MutsDos(pl+8).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);
        refere=DATA(pl+8).s(MutsDos(pl+8).Refs,2);
        srefs=nanmean(refere);
        toplot2(con,1:length(ind+8))=(DATA(pl+8).s(ind,2))-srefs;
        toplot2(toplot==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.9,.9,.9],7)
        hold on 
        text(ones(1,length(ind))*con+diferencias(i)-.003, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.01
            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('slope,s')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.120812 .12])
ylim([-.152 .13090])
ylim([-.081 .13590])
xlim([0.5 3.5])% Solo 5 cepas 
%xlim([0.5 6.5])
legend(names(orden),'Location', 'Best')

tx=text(1+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);
tx=text(2+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);
tx=text(3+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);

%% En un solo BoxPlot Tres/Dos GENES interesantes % 
figure(46)
clf
names=fieldnames(MutsDos);
orden =[13] %[5  11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1]; % platos=[3 6 2 1 4 5];
diferencias=-.375:.18:.45; %diferencias=-.375:.15:.45;
mapita=colormap(HSV);
for j=1:6;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*10,:), 'MarkerFaceColor', mapita(j*10,:) );
    hold on
end %Es para que la leyenda quede bien
 con=0;
for i = 1:length(orden)
    con=con+1;
    con2=0;
    for pl=platos;
        clear toplot toplot2
        con2=con2+1;
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        
        refere=DATA(pl).s(MutsDos(pl).Refs,2);
        srefs=nanmean(refere);
        toplot(con,1:length(ind))=(DATA(pl).s(ind,2))-srefs;
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(con2)-.055,0.001821,'patch',21,[1,1,1],mapita(con2*10,:) ,[.3 .3 .3] ,8)
        hold on
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.05
            plot( con+diferencias(con2)-.03, max(toplot(con,:))+.01, '*r')
        end
        
        refere=DATA(pl+8).s(MutsDos(pl+8).Refs,2);
        srefs=nanmean(refere);
        toplot2(con,1:length(ind))=(DATA(pl+8).s(ind,2))-srefs;
        toplot2(toplot2==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(con2)+.03,0.001821,'patch',20,[1,1,1],mapita(con2*10,:) ,[.3 .3 .3] ,8)
        hold on
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.05
            plot( con+diferencias(con2)+.03, max(toplot2(con,:))+.01, '*r')
        end
    end
    hold on
    plot( [-1 13], [0 0] )
    set(gca, 'Xtick', 1:12, 'xticklabel', [])
	%h=text( 1:12, zeros(12,1)-.151,names(orden) );
    %set( h, 'rotation', 270 )
    set(gca, 'XtickLabel',names(orden) );
    ylabel('slope,s')
end
for i=1:13
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.75 .75 .75] )
end
ylim([-.12 .065])
xlim([0.5 1.5])
legend(strains(platos),'Location', 'Best')
title('por gen')

%% En un solo BoxPlot Tres/Dos cepas interesantes
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(134)
clf
names=fieldnames(MutsDos);
orden = [5 9 8 2 1 6 13];%[5 11  12 4 3 9 8 7 2 1 10 6 13];%
platos=[1 2 5];%[3 2 5]Casos extremos de posible sesgo%[2 1]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.45:(.9/(length(orden)-1)):.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=HSV;
DATA=TBgDataAll;
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        refere=DATA(pl).s(MutsDos(pl).Refs,2);
        srefs=nanmean(refere);
        toplot(con,1:length(ind))=(DATA(pl).s(ind,2))-srefs;
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(i)-.015,.001,'patch',10,[1,1,1], mapita(i*8,:) ,[.9,.9,.9],7)
        hold on
        text(ones(1,length(ind))*con+diferencias(i)-.031, toplot(con,:), {ind} )

        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.01
            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end

        ind = MutsDos(pl+8).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);
        refere=DATA(pl+8).s(MutsDos(pl+8).Refs,2);
        srefs=nanmean(refere);
        toplot2(con,1:length(ind+8))=(DATA(pl+8).s(ind,2))-srefs;
        toplot2(toplot==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*8,:) ,[.9,.9,.9],7)
        hold on 
        text(ones(1,length(ind))*con+diferencias(i)-.003, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.01
            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('slope,s')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.120812 .12])
ylim([-.152 .13090])
ylim([-.0671 .12590])
xlim([0.5 3.5])% Solo 5 cepas 
%xlim([0.5 6.5])
%legend(names(orden),'Location', 'Best')

tx=text(1+diferencias, ones(length(diferencias),1)*-.0576, names(orden) );
set(tx, 'rotation', 270);
tx=text(2+diferencias, ones(length(diferencias),1)*-.0576, names(orden) );
set(tx, 'rotation', 270);
tx=text(3+diferencias, ones(length(diferencias),1)*-.0576, names(orden) );
set(tx, 'rotation', 270);

%%
load decayRatePLS123413y14.mat
%% Los efectos (s) correlacionan con el decayRate?
RefDecays=[];
MedRefDecays=[];
for pl=[1 2 3 4 13 14]
    w = MutsDos(pl).Refs;
    RefDecays = [RefDecays nanmean(decayRate(pl).r(1,w))]
    MedRefDecays = [MedRefDecays nanmedian(decayRate(pl).r(1,w))]
end
figure()
con=0;
for pl=[1 2 3 4 13 14]
    clear toplot    
    for i = 1:length(names)-2
        ind = MutsDos(pl).(str2mat(names(i)));
        toplot(i,1:length(ind))=DATA(pl).s(ind,2); 
    end
    con=con+1;
    %subplot(3,2,con)
    toplot(toplot==0)=NaN;
    toplot=toplot-nanmean(DATA(pl).s(MutsDos(pl).Refs,2));
    toplot=nanmedian(toplot')
    stdtoplot=std(toplot(1:12))
    plot(RefDecays(con),toplot,'.r')
	toplot=mean(toplot(1:12))
	hold on
    plot(RefDecays(con),toplot,'or')
    errorbar(RefDecays(con),toplot, stdtoplot )
    xlabel('Average Decay rate of WT')
    ylabel('Relative survival, s, of mutants')
    ylim([-.1 .1])
end 





%%
metodo = 'spearman';
%platos=[1:5 9:13];
platos=[1 9 2 10 3 11 4 12 6 14]; %acomodados por par de replica
indgenes=1:12;
genes=campos(indgenes);
mapita=[1 1 0; .9 .9 0; .8 .8 0; .7 .7 0; .6 .6 0; .5 .5 0; .4 .4 0; .3 .3 0; .15 .15 0; 0 0 0;0 .15 .15; 0 .3 .3 ; 0 .4 .4;0 .5 .5;0 .6 .6;0 .7 .7;0 .8 .8;0 .9 .9 ;0 1 1]; %amarillo/azul

aaa=clustergram(NormsMean(indgenes,platos),'RowLabels',genes,'ColumnLabels', strains(platos),'Standardize', 3, 'rowPdist', metodo,'ColumnPdist',metodo,'ColorMap', mapita  )
%'Standardize', 1, 'rowPdist', metodo,'ColumnPdist',metodo, 'Linkage','weighted', 'ColorMap', mapita);
%%
platos =[1 9 2 10 3 11 4 12 6 14];
orden=[1 2 6 10 8 9 5 11 3 4 7 12];
mapita=1-[1 1 0; .9 .9 0; .8 .8 0; .7 .7 0; .6 .6 0; .5 .5 0; .4 .4 0; .3 .3 0; .15 .15 0; 0 0 0;0 .15 .15; 0 .3 .3 ; 0 .4 .4;0 .5 .5;0 .6 .6;0 .7 .7;0 .8 .8;0 .9 .9 ;0 1 1]; %amarillo/azul
HeatMap(sMean(platos, orden), 'RowLabels', names(orden), 'ColumnLabels', strains(platos), 'ColorMap', mapita);

%% esto falta arreglar
clear NSBgDataAll temp
toremove = [Muts.RFPOnly Muts.YFPOnly];
toremove=[1 7 8];
pozos=1:20%96;
refs=0;%Muts.Refs;
mapita=JET;
%ExpBgDataAll = ExtractExponentialPoints(ODBgDataAll, 1:16, 1)
ExpBgDataAll = calculaTiempos(ExpBgDataAll, 1:16, odTh)
for platos=1:16
    col=mapita(platos*4,:)
    temp = CalcNSModel(ODBgDataAll,platos,pozos,toremove,refs,col)
    NSBgDataAll(platos)=temp(platos);
    platos
    pause
end


%%  cdfplots .s por mutante

%% Calcuar CleanData.sMean CleanData.sMedian y CleanData.std sMean y sMedian
CleanData=ODBgDataAll;
names=fieldnames(MutsDos);
%names=fieldnames(Muts);
clear sMean
for pl=1:6
    clear toplot
    for i = 1:length(names)
        %ind = Muts.(str2mat(names(i)));
        ind = MutsDos(pl).(str2mat(names(i)));
        %toplot(i,1:length(ind))=CleanData(pl).s(ind,2); 
        toplot(i,1:length(ind))= nanmean([CleanData(pl+8).s(ind,2)'; CleanData(pl).s(ind,2)']);
    end
    toplot(toplot==0)=NaN;
    toplot=toplot-nanmean(CleanData(pl).s(Muts.Refs,2));
    sMean(pl,1:15)=nanmean(toplot');
    sMedian(pl,1:15)=nanmedian(toplot');
end

%% barras por cepa
CleanData=ODBgDataAll;
names=fieldnames(MutsDos)
for i=1:8
con=0
figure()
clf
par=[1 9]+i
for pl=par % 1:16
    clear toplot    
    for i = 1:length(names)-2
        ind = MutsDos(pl).(str2mat(names(i)));
        toplot(i,1:length(ind))=CleanData(pl).s(ind,2); 
        
        
    end
    con=con+1;
    subplot(2,2,con)
    toplot(toplot==0)=NaN;
    toplot=toplot-nanmean(CleanData(pl).s(MutsDos(pl).Refs,2));
    bar(nanmedian(toplot'), 'g')
	hold on
    plot(toplot,'.r')
    title(strcat('OD-',strains(pl), ' PL-', num2str(pl)))
    ylim([-.1 .12])
end 
CleanData=TBgDataAll;

for pl=par % 1:16
    clear toplot    
    for i = 1:length(names)-2
        ind = MutsDos(pl).(str2mat(names(i)));
        toplot(i,1:length(ind))=CleanData(pl).s(ind,2); 
        
        
    end
    con=con+1;
    subplot(2,2,con)
    toplot(toplot==0)=NaN;
    toplot=toplot-nanmean(CleanData(pl).s(MutsDos(pl).Refs,2));
    bar(nanmedian(toplot'), 'g')
	hold on
    plot(toplot,'.r')
    title(strcat('T-',strains(pl), ' PL-', num2str(pl)))
    ylim([-.1 .12])
end 
end
%% AHORA POR GEN
figure(143)
CleanData=TBgDataAll;
clear pKS pW
clf
orden=[1 2 10 6 8 9 5 7 3 4 11 12];
%names=fieldnames(MutsDos);
names=fieldnames(MutsDos);
for i = 1:12%length(names)
    con=0;
    clear toplot 
	%ind = Muts.(str2mat(names(orden(i))));
    for pl=[1 9  2 10 3 11 4 12 5 13 6 14]%[3 2 6 4 1 5]%[1:6 9:14]
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        con=con+1;
        %toplot(pl,1:length(ind))=CleanData(pl).s(ind,2)-nanmedian(CleanData(pl).s(Muts.Refs,2)); 
        toplot(con,1:length(ind))=CleanData(pl).s(ind,2)-nanmedian(CleanData(pl).s(MutsDos(pl).Refs,2));
        referencias = CleanData(pl).s(MutsDos(pl).Refs,2)-nanmedian(CleanData(pl).s(MutsDos(pl).Refs,2));
        subplot(3,4,i)
        if length(find(~isnan(toplot(con,1:length(ind))))) > 0
            [h pKS(pl,i)] = kstest2(toplot(con,1:length(ind)), referencias);
            %text(con, .1, num2str(floor(p*10000)/10000) )
            hold on
            [pW(pl,i) H] = ranksum(toplot(con,1:length(ind)), referencias);
        %text(con, -.1, num2str(floor(P*10000)/10000) )
        end
    end
    toplot(toplot==0)=NaN;
    bar(nanmedian(toplot'),'g')
	hold on
    plot(toplot,'.r')
    title(strcat(names(orden(i)),'-Tmedian' ))
    ylim( [-.15 .15] )
    xlim( [0 13] )
end

%% %usar pW y pKS para solo poner las estadisticamente significativas

for i=1:12
    for pl=[1:6]%%9:14
        if pKS(i,pl)<.05 || pKS(i,pl+6)<.05
            sMeanBinKS(pl,i)=sMean(pl,i);
            %litsMeanBinKS(pl,i)=litsMean(pl,i);
        else
            sMeanBinKS(pl,i)=sMean(pl,i)*0.0;
            %litsMeanBinKS(pl,i)=litsMean(pl,i)*0.01;
        end
        if pW(i,pl)<.05 || pW(i,pl+6)<.05
            sMeanBinW(pl,i)=sMean(pl,i);
            %litsMeanBinW(pl,i)=litsMean(pl,i);
        else
            sMeanBinW(pl,i)=sMean(pl,i)*0.0;
            %litsMeanBinW(pl,i)=litsMean(pl,i)*0.01;
        end
    end
end

for i=1:12
    for pl=[1:6]%%9:14
        if pKS(i,pl)<.05 || pKS(i,pl+6)<.05 || pW(i,pl)<.05 || pW(i,pl+6)<.05
            sMeanBinKSW(pl,i)=sMean(pl,i);
            %litsMeanBinKS(pl,i)=litsMean(pl,i);
        else
            sMeanBinKSW(pl,i)=sMean(pl,i)*0.0;
            %litsMeanBinKS(pl,i)=litsMean(pl,i)*0.01;
        end
    end
end

%%
platos =[3 6 2 1 4 5];
mapita=[1 1 0; .9 .9 0; .8 .8 0; .7 .7 0; .6 .6 0; .5 .5 0; .4 .4 0; .3 .3 0; .15 .15 0; 0 0 0;0 .15 .15; 0 .3 .3 ; 0 .4 .4;0 .5 .5;0 .6 .6;0 .7 .7;0 .8 .8;0 .9 .9 ;0 1 1]; %amarillo/azul
HeatMap(sMean(platos, orden)', 'RowLabels', names(orden), 'ColumnLabels', strains(platos), 'ColorMap', mapita);
HeatMap(sMeanBinKS(platos, orden)', 'RowLabels', names(orden), 'ColumnLabels', strains(platos), 'ColorMap', mapita);
HeatMap(sMeanBinW(platos, orden)', 'RowLabels', names(orden), 'ColumnLabels', strains(platos), 'ColorMap', mapita);

%% el favorito es el spearman sMEANBn KStest solo con los pls 1:6
orden = 1:12;% [5 11  12 4 3 9 8 7 2 1 10 6];%[1 2 10 6 8 9 5 7 3 4 11 12];
platos =[3 6 2 1 4 5];
metodo = 'spearman';
mapita=[1 1 0; .9 .9 0; .8 .8 0; .7 .7 0; .6 .6 0; .5 .5 0; .4 .4 0; .3 .3 0; .15 .15 0; 0 0 0;0 .15 .15; 0 .3 .3 ; 0 .4 .4;0 .5 .5;0 .6 .6;0 .7 .7;0 .8 .8;0 .9 .9 ;0 1 1]; %amarillo/azul
mapita=[1 0 0; .9 0 0; .8 0 0; .7 0 0;.6 0 0; .5 0 0;.4 0 0; .3 0 0; .15 0  0; 0 0 0;0 0 .15; 0 0 .3 ;0 0 .4;0 0 .5;0 0 .6;0 0 .7;0 0 .8;0 0 .9 ;0 0 1]; %amarillo/azul
aaa=clustergram(sMeanBinKS(platos, orden)', 'RowLabels', names(orden), 'ColumnLabels', strains(platos), 'Standardize', 1, 'rowPdist', metodo,'ColumnPdist',metodo, 'Linkage','weighted', 'ColorMap', mapita);
addTitle(aaa, strcat(metodo, '-sMeanBinKS-OR2PLs'));




%% aplicar scripts de EGGNS
% Calculate New S,G,A con CalculaModeloNS_ScriptEGG solo los puntos Exp+T0 de cada dia
refs=Muts.Refs;
plt =1:16;
wrpl=refs;
Tiempo0=1;
ExpBgDataAll = ExtractExponentialPoints(BgDataAll, plt, 1, refs, Tiempo0);
ExpBgDataAll = calculaTiempos(ExpBgDataAll, plt, .22);
bgdataPrueba=BgDataAll2bgdataEGG(ExpBgDataAll,plt,'CFP','RFP',1);
%bgdataPrueba=ExpBgDataAll;
%BgDataSinFondo = restarFondoFPs(bgdataPrueba, [11 13],89,96)
%plt=[11 13];

%%
datExtExponential=1;
medicionesminimas=30;
extraPlRefs=0;
[bgdataAll_CMAN, data2_CMAN] = ModelASGC(bgdataPrueba,plt,wrpl,Muts.YFPOnly,Muts.RFPOnly,medicionesminimas,datExtExponential,extraPlRefs);
save 9Jul18

%%
for pl=1:8
    subplot(2,4,pl)
    %figure(pl)
    plot([-1 1],[-1 1])
    hold on    
    CorrelationScatter(TBgDataAll(pl).s(:,2),bgdataAll_CMAN(pl).S,'TiempoFijo','ModelASGC');
    p=CorrelationScatter(TBgDataAll(pl+8).s(:,2),bgdataAll_CMAN(pl+8).S,strcat('TiempoFijo PLs-',num2str([pl, pl+8]) ),'ModelASGC');
    set(p, 'color','r')
    xlim([-.15 .15])
    ylim([-.15 .15])
end


%% En un solo BoxPlot Todas juntas por cepa SNoam un box por plato
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(36)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6 13];
platos=[3 2 6 4 1 5 ]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.45:.075:.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=colormap;
DATA=bgdataAll_CMAN;%(pl).S;%TBgDataAll;
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2

        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(i)-.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        matrixmedias(i,pl) = nanmean( toplot(con,:) );
        matrixmedianas(i,pl) = nanmedian( toplot(con,:) );
        hold on
%        text(ones(1,length(ind))*con+diferencias(i)-.013, toplot(con,:), {ind} )
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere);
        if pT(orden(i),pl)<.01
%            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end
        
        ind = MutsDos(pl).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);
        toplot2(con,1:length(ind+8))=(DATA(pl+8).S(ind));
        toplot2(toplot==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        hold on        
%        text(ones(1,length(ind))*con+diferencias(i)+.013, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere);
        if pT(orden(i),pl+8)<.01
%            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
        matrixmedias2(i,pl) = nanmean( toplot2(con,:) );
        matrixmedianas2(i,pl) = nanmedian( toplot2(con,:) );

    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('slope,s')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.152 .13090])
%ylim([-.12 .085])
xlim([0.5 5.5])% Solo 5 cepas 
xlim([0.5 6.5])
legend(names(orden),'Location', 'Best')

%% En un solo BoxPlot Todas juntas por cepa SNoam promedio de los dos platos
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(65)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1 5 ]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.42:.075:.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=colormap;
DATA=bgdataAll_CMAN;%(pl).S;%TBgDataAll;
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2 toplot3
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        %notBoxPlotAA(toplot(con,:), con+diferencias(i)-.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        matrixmedias(i,pl) = nanmean( toplot(con,:) );
        matrixmedianas(i,pl) = nanmedian( toplot(con,:) );
        hold on
%        text(ones(1,length(ind))*con+diferencias(i)-.013, toplot(con,:), {ind} )
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere);
        if pT(orden(i),pl)<.01
%            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end
        
        ind = MutsDos(pl).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);
        toplot2(con,1:length(ind+8))=(DATA(pl+8).S(ind));
        toplot2(toplot==0)=NaN;
        %notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)
        hold on        
%        text(ones(1,length(ind))*con+diferencias(i)+.013, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere);
        if pT(orden(i),pl+8)<.01
%            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
        matrixmedias2(i,pl) = nanmean( toplot2(con,:) );
        matrixmedianas2(i,pl) = nanmedian( toplot2(con,:) );

        toplot3= nanmean([toplot; toplot2]);
        notBoxPlotAA(toplot(con,:), con+diferencias(i),.0025,'patch',10,[1,1,1], mapita(i*4,:) ,[.5,.5,.5],4)

    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('S, relative survivorship')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.152 .13090])
%ylim([-.12 .085])
xlim([0.5 5.5])% Solo 5 cepas 
xlim([0.5 6.5])
legend(names(orden),'Location', 'Best')

%% En un solo BoxPlot Todas juntas por gen Snoam
figure(46)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1 5]; % platos=[3 6 2 1 4 5];
diferencias=-.375:.9/length(platos):.45; %diferencias=-.375:.15:.45;
mapita=colormap(HSV);DATA=bgdataAll_CMAN;
colorpuntos=[.7,.7,.7];
saltocol= floor(64/length(platos));
for j=1:length(platos);
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*saltocol,:), 'MarkerFaceColor', mapita(j*10,:) );
    hold on
end %Es para que la leyenda quede bien
 con=0;
for i = 1:length(orden)
    con=con+1;
    con2=0;
    for pl=platos;
        clear toplot toplot2
        con2=con2+1;
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(con2)-.03,0.001,'patch',20,[1,1,1],mapita(con2*saltocol,:) ,colorpuntos,3)
        hold on
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere);
        if pT(orden(i),pl)<.05
            plot( con+diferencias(con2)-.03, max(toplot(con,:))+.01, '*r')
        end
        
        toplot2(con,1:length(ind))=(DATA(pl+8).S(ind));
        toplot2(toplot2==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(con2)+.03,0.001,'patch',20,[1,1,1],mapita(con2*saltocol,:) ,colorpuntos,3)
        hold on
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere);
        if pT(orden(i),pl+8)<.05
            plot( con+diferencias(con2)+.03, max(toplot2(con,:))+.01, '*r')
        end
    end
    hold on
    plot( [-1 13], [0 0] )
    set(gca, 'Xtick', 1:12, 'xticklabel', [])
	%h=text( 1:12, zeros(12,1)-.151,names(orden) );
    %set( h, 'rotation', 270 )
    set(gca, 'XtickLabel',names(orden) );
    ylabel('slope,s')
end
for i=1:13
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.75 .75 .75] )
end
ylim([-.172 .1175])
xlim([0.5 12.5])
legend(strains(platos),'Location', 'Best')
title('por gen')

%% En un solo BoxPlot Todas juntas por gen Snoam
figure(47)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1 5]; % platos=[3 6 2 1 4 5];
diferencias=-.375:.9/length(platos):.45; %diferencias=-.375:.15:.45;
mapita=colormap(HSV);DATA=bgdataAll_CMAN;
colorpuntos=[.7,.7,.7];
saltocol= floor(64/length(platos));
for j=1:length(platos);
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*saltocol,:), 'MarkerFaceColor', mapita(j*10,:) );
    hold on
end %Es para que la leyenda quede bien
 con=0;
for i = 1:length(orden)
    con=con+1;
    con2=0;
    for pl=platos;
        clear toplot toplot2
        con2=con2+1;
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(con2)-.03,0.001,'patch',20,[1,1,1],mapita(con2*saltocol,:) ,colorpuntos,3)
        hold on
        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere);
        if pT(orden(i),pl)<.05
            plot( con+diferencias(con2)-.03, max(toplot(con,:))+.01, '*r')
        end
        
        toplot2(con,1:length(ind))=(DATA(pl+8).S(ind));
        toplot2(toplot2==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(con2)+.03,0.001,'patch',20,[1,1,1],mapita(con2*saltocol,:) ,colorpuntos,3)
        hold on
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere);
        if pT(orden(i),pl+8)<.05
            plot( con+diferencias(con2)+.03, max(toplot2(con,:))+.01, '*r')
        end
    end
    hold on
    plot( [-1 13], [0 0] )
    set(gca, 'Xtick', 1:12, 'xticklabel', [])
	%h=text( 1:12, zeros(12,1)-.151,names(orden) );
    %set( h, 'rotation', 270 )
    set(gca, 'XtickLabel',names(orden) );
    ylabel('slope,s')
end
for i=1:13
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.75 .75 .75] )
end
ylim([-.172 .1175])
xlim([0.5 12.5])
legend(strains(platos),'Location', 'Best')
title('por gen')



%% En un solo BoxPlot Tres/Dos cepas interesantes
%nombres={'Clinical-1, Italy','Clinical-2, Italy','Wine, Chile','Laboratory
%(s288c)','Wine, France','Woodlands, USA'}
figure(34)
clf
names=fieldnames(MutsDos);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];%[5 9 8 2 1 6 13]%
platos=[1 2 5]%[3 2 5]Casos extremos de posible sesgo%[2 1]; %platos=[3 14 2 1 4 13]; %platos=[3 6 2 1 4 5];
diferencias=-.45:.075:.45;%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=colormap;
DATA=bgdataAll_CMAN;%(pl).S;%TBgDataAll;
colorpuntos=[.7,.7,.7];
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        %ind = Muts.(str2mat(names(orden(i))));
        con=con+1;

        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        notBoxPlotAA(toplot(con,:), con+diferencias(i)-.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,colorpuntos,7)
        hold on
%        text(ones(1,length(ind))*con+diferencias(i)-.031, toplot(con,:), {ind} )

        [hT( orden(i),pl ) pT(orden(i),pl)] = ttest2(toplot(con,:), refere-srefs);
        if pT(orden(i),pl)<.01
            plot( con+diferencias(i)-.013, max(toplot(con,:))+.01, '*r')
        end

        ind = MutsDos(pl+8).(str2mat(names(orden(i))));
     %   toplot2(con,1:length(ind))=CleanData(pl+8).snorm(ind);

        toplot2(con,1:length(ind+8))=(DATA(pl+8).S(ind));
        toplot2(toplot==0)=NaN;
        notBoxPlotAA(toplot2(con,:), con+diferencias(i)+.013,.001,'patch',10,[1,1,1], mapita(i*4,:) ,colorpuntos,7)
        hold on 
%        text(ones(1,length(ind))*con+diferencias(i)-.003, toplot2(con,:), {ind} )
        [hT( orden(i),pl+8 ) pT(orden(i),pl+8)] = ttest2(toplot2(con,:), refere-srefs);
        if pT(orden(i),pl+8)<.01
            plot( con+diferencias(i)+.013, max(toplot2(con,:))+.01, '*r')
        end
    end
    %(y,x,jitter,style, jitterBox, sdcolor, SEMcolor, pointColor, MarkerSize)
    hold on
    plot( [-1 7], [0 0] )
%    set(gca, 'Xtick', 1:length(strains(platos) ), 'xticklabel', [])
	%h=text
    set(gca,'Xtick',1:length(platos),'XtickLabel', strains(platos) );
%    set( h, 'rotation', 270 )
    ylabel('S, ro-rx')
end
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.120812 .12])
ylim([-.152 .13090])
ylim([-.081 .13590])
xlim([0.5 3.5])% Solo 5 cepas 
%xlim([0.5 6.5])
legend(names(orden),'Location', 'Best')

tx=text(1+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);
tx=text(2+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);
tx=text(3+diferencias, ones(length(diferencias),1)*-.0676, names(orden) );
set(tx, 'rotation', 270);


%%
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
platos=[3 2 6 4 1 5]; % platos=[3 6 2 1 4 5];
mapita=[1 1 0; .9 .9 0; .8 .8 0; .7 .7 0; .6 .6 0; .5 .5 0; .4 .4 0; .3 .3 0; .15 .15 0; 0 0 0;0 .15 .15; 0 .3 .3 ; 0 .4 .4;0 .5 .5;0 .6 .6;0 .7 .7;0 .8 .8;0 .9 .9 ;0 1 1]; %amarillo/azul
metodo='spearman';
clustergram(promediomatrixmedias(1:12,:),'RowLabels',genes(orden),'ColumnLabels', strains(platos),'Standardize', 3, 'rowPdist', metodo,'ColumnPdist',metodo,'ColorMap', mapita  )


%% Plot por cepa, colores por cepa y xticklabels con nombres de genes
%load CleanDataOaxaca2017.mat
figure(332)
clf
names=fieldnames(Muts);
orden = [5 11  12 4 3 9 8 7 2 1 10 6];
%platos=[3 14 2 1 4 13];
platos=[3 6 2 1 4 5]
diferencias=-.425:.075:.45%[-.42 -.34 -.26 -.18 -.1 -.05 .05 .1 .18 .26 .34 .42 ];
mapita=HSV;
clf
ToXTick=[];
for j=1:13;
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*4,:), 'MarkerFaceColor', mapita(j*4,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=0;
    for pl=platos;
        clear toplot toplot2
        ind = Muts.(str2mat(names(orden(i))));
        con=con+1;
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        toplot2(con,1:length(ind))=(DATA(pl+8).S(ind));
        toplot2(toplot==0)=NaN;
        ToplotMean=nanmean([toplot(con,:) toplot2(con,:)]);
        ToplotStd=nanstd([toplot(con,:) toplot2(con,:)]/sqrt(12) );
        PlotConError(con+diferencias(i),ToplotMean, ToplotMean+ToplotStd, ToplotMean-ToplotStd,mapita(con*10,:),'*',10,.01 )
        ToXTick=[ToXTick con+diferencias(i) ]
%        PlotConError(con+diferencias(i),ToplotMean, ToplotMean+ToplotStd, ToplotMean-ToplotStd,mapita(i*4,:),'o',6,.01 )
    end

end
    hold on
    plot( [-1 7], [0 0] )
    %set(gca, 'Xtick', 1:6, 'xticklabel', [])
	h=text( [.6 1.6 2.56 3.6 4.51 5.58], zeros(6,1)+.108,strains(platos), 'fontsize', 14 );
    set( h, 'rotation', 0 )
    ylabel('S, relative survival')
for i=1:7
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.5 .5 .5] )
end
ylim([-.12 .1])
xlim([0.5 6.5])
xlabels=[names(orden)];
xlabels=[xlabels;xlabels;xlabels;xlabels;xlabels;xlabels;];
xcoordinates=[0.575:.075:1.4];%-.425:.075:.45
xcoordinates=[xcoordinates xcoordinates+1.01 xcoordinates+2.01 xcoordinates+3.01 xcoordinates+4.01 xcoordinates+5.01];
%h=text( xcoordinates, zeros(length(xcoordinates),1)-.121, lower(xlabels),'fontsize',14)
%set(h, 'rotation',270)
set(gca, 'xtick',sort(ToXTick), 'xticklabel', {ones(6,1)*1:12})

%% Plot por cepa, colores por cepa y xticklabels con nombres de genes En un solo BoxPlot Todas juntas por gen Snoam
figure(47)
clf
names=fieldnames(MutsDos);
orden = [ 7 2 1 10 6];%5 11  12 4 3 9 8
platos=[3 2 6 4 1 5]; % platos=[3 6 2 1 4 5];
diferencias=-.375:.9/length(platos):.45; %diferencias=-.375:.15:.45;
mapita=colormap(HSV);DATA=bgdataAll_CMAN;
colorpuntos=[.7,.7,.7];
saltocol= floor(64/length(platos));
for j=1:length(platos);
    p=plot(-5,-5,'sb','MarkerSize',7);
    set(p,'Color', mapita(j*saltocol,:), 'MarkerFaceColor', mapita(j*10,:) );
    hold on
end %Es para que la leyenda quede bien
for i = 1:length(orden)
    con=con+1;
    con2=0;
    for pl=platos;
        clear toplot toplot2
        

        
        con2=con2+1;
        ind = MutsDos(pl).(str2mat(names(orden(i))));
        
        toplot(con,1:length(ind))=(DATA(pl).S(ind));
        toplot(toplot==0)=NaN;
        toplot2(con,1:length(ind))=(DATA(pl+8).S(ind));
        toplot2(toplot2==0)=NaN;
        
        ToplotMean=nanmean([toplot(con,:) toplot2(con,:)]);
        ToplotStd=nanstd([toplot(con,:) toplot2(con,:)]/sqrt(12) );
        PlotConError(i+diferencias(con2),ToplotMean, ToplotMean+ToplotStd, ToplotMean-ToplotStd,mapita(con2*10,:),'*',10,.01 )

    end
    hold on
    plot( [-1 10], [0 0] )
    set(gca, 'Xtick', 1:12, 'xticklabel', [])
	%h=text( 1:12, zeros(12,1)-.151,names(orden) );
    %set( h, 'rotation', 270 )
    set(gca, 'XtickLabel',names(orden) );

end
for i=1:13
    plot([i-.5 i-.5], [1.8 -.8], '--', 'Color', [.75 .75 .75] )
end
 con=0;
 ylim([-.172 .1175])
xlim([0.5 length(orden)+.5])
legend(strains(platos),'Location', 'Best')
title('por gen')
ylabel('S, relative survival')

%%
% para figuras del congreso heidelberg

wells=Muts.Tpk1;
refs=Muts.Refs;
figure(22)
clf
con=0;
for pl= [19 17]
    con=con+1;
    subplot( 1,2, con )
    x = nanmean(survival(pl).t(:,wells)') ;
    y = nanmean(survival(pl).s(:,wells)') ;
    stdy = nanstd(survival(pl).s(:,wells)')./sqrt(length(wells)) ;
    plotconError(x, y, y-stdy, y+stdy,mapita(10*(con+2),:),'o-', 1, .001 );
    hold on
    x = nanmean(survival(pl).t(:,refs)') ;
    y = nanmean(survival(pl).s(:,refs)') ;
    stdy = nanstd(survival(pl).s(:,refs)')./sqrt(length(refs)) ;
    plotconError(x, y, y-stdy, y+stdy,'k','o--', 1, .001 );    
    title( num2str(pl) )
    
    xlim( [-.5 12] )
    ylim( [10 105] )
    set(gca, 'ytick', 25:25:100, 'xtick', 03:12)
end






