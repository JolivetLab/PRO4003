
% Initiate directory for saving data.
thisDirectory   = fileparts(mfilename('fullpath'));
saveDirectory   = fullfile(thisDirectory,'Cullen2018_R1_Figure_S8');
if ~isdir(saveDirectory)
    mkdir(saveDirectory)
end

% Regenerate MAT files for the channels.
activeChannel   = McIntyre2002SlowK;        save('SavedParameters/ActiveChannels/McIntyre2002SlowK.mat','activeChannel');
activeChannel   = McIntyre2002FastNa;       save('SavedParameters/ActiveChannels/McIntyre2002FastNa.mat','activeChannel');
activeChannel   = McIntyre2002PersistentNa; save('SavedParameters/ActiveChannels/McIntyre2002PersistentNa.mat','activeChannel');
clear activeChannel

% Initiate temperature.
temp            = [21 37];

% Figure.
f               = figure;

% Already adding panels for insets.
ax1 = axes('Position',[0.9 0.25 0.05 0.05]);
ax2 = axes('Position',[0.9 0.20 0.05 0.05]);
ax3 = axes('Position',[0.9 0.15 0.05 0.05]);
ax4 = axes('Position',[0.9 0.10 0.05 0.05]);

% Run the model and calculate CV.
for k = 1:2
    
    % Produce parameters for default cortex model.
    clear par;
    par = Cullen2018CortexAxon();
    
    % Set temperature.
    par.sim.temp = temp(k);
    
    %% Run controls for insets.
    for psw = [0 20]
      	par.sim.dt.value                    = 2;
        par.sim.tmax.value                  = 800;
        par.myel.geo.peri.value.ref         = psw;
        par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                                 = CalculateNumberOfMyelinLamellae(par, 'max');
        par.node.geo.length.value.ref       = 1.492*exp(-psw/4.918)+0.4;
        par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
        par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
        par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
        par =                                 CalculateLeakConductance(par);
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, fullfile(saveDirectory, ['Cullen2018Cortex_CTRL_psw_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
        if k == 1  && psw == 0
            plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-c','Parent',ax3);
        elseif k == 1  && psw == 20
            plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-c','Parent',ax4);
        elseif k == 2  && psw == 0
            plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-r','Parent',ax1); 
        else
            plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-r','Parent',ax2);
        end
    end
    refresh;
   
    % Reset model and temperature.
    clear par;
    par                     = Cullen2018CortexAxon();
    par.sim.temp            = temp(k);
    par.sim.tmax.value      = 15;

    %% Run all simulations varying periaxonal space width.
    j = 1;
    for psw = [0:0.2:1.6 2:6 6.477 7:8 8.487 9 10:2:14 15 20]
        par.myel.geo.peri.value.ref         = psw;
        par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                                 = CalculateNumberOfMyelinLamellae(par, 'max'); 
        par.node.geo.length.value.ref       = 1.492*exp(-psw/4.918)+0.4;
        par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
        par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
        par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
        par =                                 CalculateLeakConductance(par);
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, fullfile(saveDirectory, ['Cullen2018Cortex_psw_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
        velocity.psw(j,1)           = psw;
        if k == 1
            velocity.psw(j,2)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        else
            velocity.psw(j,3)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
        j                           = j+1;
    end
    refresh;
    
end

% Finish writing.
dlmwrite([saveDirectory '/psws.txt'],velocity.psw);

% Finish plot.
subplot(131);
plot(velocity.psw(:,1),velocity.psw(:,2),'-c'), hold on
plot(velocity.psw(:,1),velocity.psw(:,3),'-r');

% Add basics.
subplot(131); axis([0 20 0 4.5]), xticks(0:5:20), xticklabels({'0','5','10','15','20'}), yticks(0:1:4), yticklabels({'0','1','2','3','4'}), xlabel('Periaxonal space width (nm)'), ylabel('Conduction velocity (m/s)');
hold on, plot([6.477 6.477],[0 4.5],'.-k'); set(gca,'Box','off');

% Add last graph.
Chosen_d = 1e-02;
subplot(132);
plot(velocity.psw(:,1),1e03*Chosen_d./velocity.psw(:,2),'-c'); hold on;
plot(velocity.psw(:,1),1e03*Chosen_d./velocity.psw(:,3),'-r');
axis([0 20 0 15]);
xticks(0:5:20), xticklabels({'0','5','10','15','20'}), xlabel('Periaxonal space width (nm)');
yticks(0:5:15), yticklabels({'0','5','10','15'}), ylabel('Conduction delay over 1 cm (ms)');
plot([6.477 6.477],[0 15],'.-k');
text(10,13.5,['\Deltat_{0 - 20} = ~' num2str(1e03*Chosen_d./velocity.psw(end,2)-1e03*Chosen_d./velocity.psw(1,2),'%1.0f') ' ms/cm'],'Color','c');
text(10,3, ['\Deltat_{0 - 20} = ~' num2str(1e03*Chosen_d./velocity.psw(end,3)-1e03*Chosen_d./velocity.psw(1,3),'%1.0f') ' ms/cm'],'Color','r');

% Finish insets.
uistack(ax1,'top'); uistack(ax2,'top'); uistack(ax3,'top'); uistack(ax4,'top');
set(f,'CurrentAxes',ax1); axis([0 200 -90 -40]); set(gca,'Box','off','Position',[0.28 0.27 0.05 0.05]); xticklabels({});
set(f,'CurrentAxes',ax2); axis([0 200 -90 -40]); set(gca,'Box','off','Position',[0.34 0.27 0.05 0.05]); xticklabels({}); yticklabels({});
set(f,'CurrentAxes',ax3); axis([0 200 -90 -40]); set(gca,'Box','off','Position',[0.28 0.20 0.05 0.05]);
set(f,'CurrentAxes',ax4); axis([0 200 -90 -40]); set(gca,'Box','off','Position',[0.34 0.20 0.05 0.05]); yticklabels({});
subplot(132); set(gca,'Box','off','Position',get(gca,'Position')+[0.05 0 0 0]);

% Printing as JPEG.
set(gcf,...
    'PaperUnits','inches',...
    'PaperSize',[15 15],...
    'PaperPosition',[0.5 0.5 14 14],...
    'Renderer','Painters');
print(gcf,[saveDirectory '/Cullen2018_R1_Figure_S8.jpg'],'-djpeg','-r300');
print(gcf,[saveDirectory '/Cullen2018_R1_Figure_S8.pdf'],'-dpdf');
