%%
% The code used to generate all 7 figures are in one .m file, the code
% ...should be run in sections to display all figures in different windows.
% The figures may slightly different in higher version of Matlab
%% Calculation tables of net CO2 emissions
clc,clear,close all
load ('Net_CO2_tables');

%% Q1.1
clc,clear,close all
x = linspace(0,1); % define the scale of x axis(efficiency)(0 to 100%)
y1 = 0.276485044*ones(size(x)); % net CO2 emission of landfill
y2 = 0.45562 - 1.19969109*x; % net CO2 emission of 
plot(x, y2,':k')
hold on
plot(x, y1,'--');
[xi,yi] = polyxpoly(x,y1,x,y2); % display the intersection point of two line
xlabel('\eta of EfW')
ylabel('Net tCO2 or CO2e emitted') 
legend('EfW','Landfill')
grid on
disp([xi,yi]);
%% Q1.2/1.3
clc, clear, close all
n = 1;
while n <= 1000; % set the sample number
    C = rand(1,11); % randomly generate 11 categories of waste as a vector
    C = C/sum(C); % define the total value of 100%
    factor = 0.373; % define grid intensity
    Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45]; % vector of biogenic C proportion
    Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0]; % vector of fossil C proportion
    d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13]; % vector of decomposable C proportion
    c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1]; % vector of calorific values
    % EfW part
    m_fos_c = C.*Fos; % mass of fossil carbon for each category 
    m_bio_c = C.*Bio; % mass of biogenic carbon for each category
    fos_CO2 = m_fos_c*(44/12); % mass of fossil CO2 for each category
    bio_CO2 = m_bio_c*(44/12); % mass of fossil CO2 for each category
    E = c_v.*C; % calorific value for each category 
    efw_CCGT_offset = sum(factor*E); % total CO2 offset of CCGT
    Net_efw = sum(fos_CO2) - efw_CCGT_offset; % net CO2 emission of EfW
    % Landfill part
    CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28; % CO2e_CH4 of landfill
    landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor; %CCGT CO2 offset of landfill
    Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset); % net CO2 emission of landfill
    eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset; %Minimum EfW efficiency
    x(n) = m_bio_c/(m_bio_c+m_fos_c); % define x axis 
    n = n+1;
end
scatter(x(1:n-1),eff(1:n-1),'b'); % plot the result point of each sample
hold on % preserve all points in one graph
line([0 1],[0 0],'Color','Black') % draw a line y=0
hold on
p = polyfit(x,eff,1); % coefficients of fit line 
v = polyval(p,x); % draw fit line on the graph
plot(x,v,'k');
xlabel('Biogenic carbon proportion');
ylabel('Minimum EfW efficiency \eta')
hold all % diplay all lines and points
grid on
%% 2.1 carbon intensity (random generated)
clc,clear, close all
for i = 1:10 % define grid intensity sample number as 10
    factor = rand; % randomly generate the factor (0 to 1)
    % waste components generation same as above, but as nested loop
    n = 1;
    while n <= 100;
        C = rand(1,11);
        C = C/sum(C);
        Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
        Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
        d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
        c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
        m_fos_c = C.*Fos;
        m_bio_c = C.*Bio;
        fos_CO2 = m_fos_c*(44/12);
        bio_CO2 = m_bio_c*(44/12);
        E = c_v.*C;
        efw_CCGT_offset = sum(factor*E);
        Net_efw = sum(fos_CO2) - efw_CCGT_offset;
        
        CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
        landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor;
        Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
        eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
        x(n) = m_bio_c/(m_bio_c+m_fos_c);
        n = n+1;
    end
    %scatter(x(1:n-1),eff(1:n-1),'b');
    %hold on
    p = polyfit(x,eff,1);
    v = polyval(p,x);
    txt = ['carbon intensity:',num2str(factor)];% generate legends in different colors and names
    plot(x,v,'display',txt); % generate legends in different colors and names
    xlabel('Biogenic carbon proportion');
    ylabel('Minimum EfW efficiency \eta')
    hold all
    grid on
end
legend show
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)
%% carbon intensity (yearly scenario)
clc,clear, close all
% IN year 2008
factor_08 = 0.353*0.373+0.341*0.9+0.061*0.6
n = 1;
while n <= 100;
    C = rand(1,11);
    C = C/sum(C);
    Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
    Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
    d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
    c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
    m_fos_c = C.*Fos;
    m_bio_c = C.*Bio;
    fos_CO2 = m_fos_c*(44/12);
    bio_CO2 = m_bio_c*(44/12);
    E = c_v.*C;
    efw_CCGT_offset = sum(factor_08*E);
    Net_efw = sum(fos_CO2) - efw_CCGT_offset;
    
    CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
    landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor_08;
    Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
    eff1(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
    x1(n) = m_bio_c/(m_bio_c+m_fos_c);
    n = n+1;
end
%scatter(x(1:n-1),eff(1:n-1),'b');
%hold on
p1 = polyfit(x1,eff1,1);
v1 = polyval(p1,x1);
plot(x1,v1); % generate legends in different colors and names
xlabel('Biogenic carbon proportion');
ylabel('Minimum EfW efficiency \eta')
hold all
grid on
% IN year 2013
factor_13 = 0.386*0.373+0.215*0.9+0.02*0.6
n = 1;
while n <= 100;
    C = rand(1,11);
    C = C/sum(C);
    Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
    Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
    d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
    c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
    m_fos_c = C.*Fos;
    m_bio_c = C.*Bio;
    fos_CO2 = m_fos_c*(44/12);
    bio_CO2 = m_bio_c*(44/12);
    E = c_v.*C;
    efw_CCGT_offset = sum(factor_13*E);
    Net_efw = sum(fos_CO2) - efw_CCGT_offset;
    
    CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
    landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor_13;
    Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
    eff2(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
    x2(n) = m_bio_c/(m_bio_c+m_fos_c);
    n = n+1;
end
%scatter(x(1:n-1),eff(1:n-1),'b');
%hold on
p2 = polyfit(x2,eff2,1);
v2 = polyval(p2,x2);
plot(x2,v2); % generate legends in different colors and names
xlabel('Biogenic carbon proportion');
ylabel('Minimum EfW efficiency \eta')
hold all
grid on
% IN year 2018
factor_18 = 0.335*0.373+0.11*0.9+0.004*0.6
n = 1;
while n <= 100;
    C = rand(1,11);
    C = C/sum(C);
    Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
    Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
    d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
    c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
    m_fos_c = C.*Fos;
    m_bio_c = C.*Bio;
    fos_CO2 = m_fos_c*(44/12);
    bio_CO2 = m_bio_c*(44/12);
    E = c_v.*C;
    efw_CCGT_offset = sum(factor_18*E);
    Net_efw = sum(fos_CO2) - efw_CCGT_offset;
    
    CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
    landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor_18;
    Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
    eff3(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
    x3(n) = m_bio_c/(m_bio_c+m_fos_c);
    n = n+1;
end
%scatter(x(1:n-1),eff(1:n-1),'b');
%hold on
p3 = polyfit(x3,eff3,1);
v3 = polyval(p3,x3);
plot(x3,v3); % generate legends in different colors and names
xlabel('Biogenic carbon proportion');
ylabel('Minimum EfW efficiency \eta')
hold all
grid on
legend('2008','2013','2018');
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)
%% 2.2 CH4 capture proportion
clc,clear, close all
for i = 1:10 % incremental steps of CH4 capture proportion change for 10 times
    c_rate = i*0.1; %scale from 0.1 to 1
    % the original capture rate was replaced by c_rate in nested loop
    n = 1;
    % waste components generation same as above, but as nested loop
    while n <= 50;
        C = rand(1,11);
        C = C/sum(C);
        factor = 0.373;
        Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
        Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
        d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
        c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
        m_fos_c = C.*Fos;
        m_bio_c = C.*Bio;
        fos_CO2 = m_fos_c*(44/12);
        bio_CO2 = m_bio_c*(44/12);
        E = c_v.*C;
        efw_CCGT_offset = sum(factor*E);
        Net_efw = sum(fos_CO2) - efw_CCGT_offset;
        
        CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
        landfill_CCGT_offset = C.*d_c*0.5*(16/12)*c_rate*0.5*2.84*factor;
        Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
        eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
        x(n) = m_bio_c/(m_bio_c+m_fos_c);
        n = n+1;
    end
    %scatter(x(1:n-1),eff(1:n-1),'b');
    %hold on
    p = polyfit(x,eff,1);
    v = polyval(p,x);
    txt = ['captured:',num2str(c_rate)];
    plot(x,v,'display',txt);
    xlabel('Biogenic carbon proportion');
    ylabel('Minimum EfW efficiency \eta')
    hold all
    grid on
end
legend show
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)
%% 2.2 Proportion of decomposable carbon that becomes methane
clc,clear, close all
for i = 0:10 % incremental steps of decomposes to CH4 proportion for 10 times
    factor = 0.373;
    d_CH4 = 0.1*i; % scale from 0.1 to 1
    n = 1;
    while n <= 100;
        C = rand(1,11);
        C = C/sum(C);
        Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
        Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
        d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
        c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
        m_fos_c = C.*Fos;
        m_bio_c = C.*Bio;
        fos_CO2 = m_fos_c*(44/12);
        bio_CO2 = m_bio_c*(44/12);
        E = c_v.*C;
        efw_CCGT_offset = sum(factor*E);
        Net_efw = sum(fos_CO2) - efw_CCGT_offset;
        
        CO2e_CH4 = C.*d_c*d_CH4*(16/12)*0.25*0.9*28;
        landfill_CCGT_offset = C.*d_c*d_CH4*(16/12)*0.75*0.5*2.84*factor;
        Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
        eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
        x(n) = m_bio_c/(m_bio_c+m_fos_c);
        n = n+1;
    end
    %scatter(x(1:n-1),eff(1:n-1),'b');
    %hold on
    p = polyfit(x,eff,1);
    v = polyval(p,x);
    txt = ['decomp to CH4:',num2str(d_CH4)];
    plot(x,v,'display',txt);
    xlabel('Biogenic carbon proportion');
    ylabel('Minimum EfW efficiency \eta')
    hold all
    grid on
end
legend show
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)

%% proportion of CH4 oxidized
clc,clear, close all
for i = 0:10 % incremental steps of oxidized CH4 proportion
    factor = 0.373;
    o_CH4 = 0.1*i; % scale from 10% to 100%
    n = 1;
    while n <= 100;
        C = rand(1,11);
        C = C/sum(C);
        Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
        Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
        d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
        c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
        m_fos_c = C.*Fos;
        m_bio_c = C.*Bio;
        fos_CO2 = m_fos_c*(44/12);
        bio_CO2 = m_bio_c*(44/12);
        E = c_v.*C;
        efw_CCGT_offset = sum(factor*E);
        Net_efw = sum(fos_CO2) - efw_CCGT_offset;
        
        CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*(1-o_CH4)*28;
        landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor;
        Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
        eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
        x(n) = m_bio_c/(m_bio_c+m_fos_c);
        n = n+1;
    end
    %scatter(x(1:n-1),eff(1:n-1),'b');
    %hold on
    p = polyfit(x,eff,1);
    v = polyval(p,x);
    txt = ['oxidized CH4:',num2str(o_CH4)];
    plot(x,v,'display',txt);
    xlabel('Biogenic carbon proportion');
    ylabel('Minimum EfW efficiency \eta')
    hold all
    grid on
end
legend show
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)
%% Calorific Value (CV) of waste
clc,clear, close all
for i = 0:10 % incremental steps of total CV
    factor = 0.373;
    times = 0.3*i; % times of original total CV from 0.3 to 3
    n = 1;
    while n <= 100;
        C = rand(1,11);
        C = C/sum(C);
        Bio = [0.32, 0, 0.2, 0.19, 0.04, 0.14, 0.17, 0.07, 0, 0, 0.45];
        Fos = [0, 0.52, 0.2, 0.19, 0.04, 0, 0, 0, 0, 0, 0];
        d_c = [0.16, 0, 0.07, 0.09, 0, 0.09, 0.09, 0, 0, 0, 0.13];
        c_v = [3.5, 7.05, 4.44, 4.33, 0.78, 1.47, 1.85, 1.3, 0.43, 0, 5.1];
        m_fos_c = C.*Fos;
        m_bio_c = C.*Bio;
        fos_CO2 = m_fos_c*(44/12);
        bio_CO2 = m_bio_c*(44/12);
        E = times*c_v.*C;
        efw_CCGT_offset = sum(factor*E);
        Net_efw = sum(fos_CO2) - efw_CCGT_offset;
        
        CO2e_CH4 = C.*d_c*0.5*(16/12)*0.25*0.9*28;
        landfill_CCGT_offset = C.*d_c*0.5*(16/12)*0.75*0.5*2.84*factor;
        Net_landfill = sum(CO2e_CH4 - landfill_CCGT_offset);
        eff(n) = (sum(fos_CO2) - Net_landfill)/efw_CCGT_offset;
        x(n) = m_bio_c/(m_bio_c+m_fos_c);
        n = n+1;
    end
    %scatter(x(1:n-1),eff(1:n-1),'b');
    %hold on
    p = polyfit(x,eff,1);
    v = polyval(p,x);
    txt = ['CV times:',num2str(times)];
    plot(x,v,'display',txt);
    xlabel('Biogenic carbon proportion');
    ylabel('Minimum EfW efficiency \eta')
    hold all
    grid on
end
legend show
y0 = '\eta=0';
line([0 1],[0 0],'Color','Black','display',y0)