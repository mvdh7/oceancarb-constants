% Define input conditions
dic = (1600:3200)';
ta = 2250;
t = 25;
s = 35;

% Evaluate buffer factors
[gDIC, bDIC, oDIC, gALK, bALK, oALK] = ESM10_CO2SYS( ...
    ta, dic, 1, 2, s, t, t, 0, 0, 0, 0, 3, 10, 3);

% Plot results like ESM10 Fig. 2
figure(1); clf; hold on

plot(dic/1000, gDIC/1000, 'linewidth',1, 'color',[0 0 0.5])
plot(dic/1000, bDIC/1000, 'linewidth',1, 'color',[0.5 0 0.5])
plot(dic/1000,-oDIC/1000, 'linewidth',1, 'color',[0.8 0 0])
plot(dic/1000,-gALK/1000, 'linewidth',1, 'color',[0.5 0 0.5])
plot(dic/1000,-bALK/1000, 'linewidth',1, 'color',[0 0.4 0])
plot(dic/1000, oALK/1000, 'linewidth',1, 'color',[0 0.6 0.6])

xlim([1.6 3.2])
ylim([0.1 1])
set(gca, 'box','on', 'xtick',1.6:0.2:3.2, 'ytick',0.1:0.1:1)
xlabel('DIC (mM)')
ylabel('buffer factor (mM)')

legend('\gamma_{DIC}','-\gamma_{Alk}','\beta_{DIC}','-\beta_{Alk}', ...
    '-\omega_{DIC}','\omega_{Alk}', 'location','nw')
