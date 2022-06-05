
% Initiate directory for saving data.
thisDirectory   = fileparts(mfilename('fullpath'));
saveDirectory   = fullfile(thisDirectory,'Cullen2018_R0_Convergence');
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

% Run the model and calculate CV.
for k = 1:2
    
    % Index.
    ns = 1;
    
    for nseg = [1 2 3 4 5 6 7 10 15 20 30 40 50 52]
        
        
        % Produce parameters for default cortex model.
        clear par;
        par = Cullen2018CortexAxon();
        
        % Set temperature.
        par.sim.temp                        = temp(k);
        
        % Adjust simulation time.
        par.sim.tmax.value                  = 15;
        
        % Adjust number of segments.
        par.geo.nintseg                     = nseg;
        par.intn.seg.geo.length.value.ref   = par.intn.geo.length.value.ref / par.geo.nintseg;
        par.intn.seg.geo.length.value.vec   = repmat(par.intn.geo.length.value.vec / par.geo.nintseg, 1, par.geo.nintseg);
        par.intn.seg.geo.length.units       = {1, 'um', 1};
        
        % Internode segment diameter (=internode diameter).
        par.intn.seg.geo.diam.value.ref     = par.intn.geo.diam.value.ref;
        par.intn.seg.geo.diam.value.vec     = repmat(par.intn.geo.diam.value.vec, 1, par.geo.nintseg);
        par.intn.seg.geo.diam.units         = {1, 'um', 1};
        
        % Periaxonal space width.
        par.myel.geo.peri.value.ref         = 6.477;
        par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.peri.units             = {1, 'nm', 1};
        
        % Myelin wrap periodicity.
        par.myel.geo.period.value           = 16.2873;
        par.myel.geo.period.units           = {1, 'nm', 1};
        
        % g-ratio (internode axon diameter to internode outer diameter ratio)
        par.myel.geo.gratio.value.ref       = 0.724;
        par.myel.geo.gratio.value.vec_ref   = par.myel.geo.gratio.value.ref * ones(par.geo.nintn, par.geo.nintseg);
        
        
        
        %% Run sham simulations.
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, fullfile(saveDirectory, ['Cullen2018Cortex_sham_' num2str(par.sim.temp) 'C.mat']));
        if k == 1
            velocity.lotemp.nseg(ns,1) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        else
            velocity.hitemp.nseg(ns,1) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
        dlmwrite([saveDirectory '/time_vector_sham_' num2str(par.sim.temp) 'C.txt'],TIME_VECTOR);
        dlmwrite([saveDirectory '/membrane_potential_sham_' num2str(par.sim.temp) 'C.txt'],MEMBRANE_POTENTIAL);
        refresh;
        
        %% Run simulations for extreme periaxonal space widths.
        j = 2;
        for psw = [0 20]
            par.myel.geo.peri.value.ref = psw;
            par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
            par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
            par                         = CalculateNumberOfMyelinLamellae(par, 'max');
            [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par, fullfile(saveDirectory, ['Cullen2018Cortex_psw_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
            if k == 1
                velocity.lotemp.nseg(ns,j)  = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
            else
                velocity.hitemp.nseg(ns,j)	= velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
            end
            j = j+1;
        end
        
        %% Update counter.
        ns = ns+1;
        
    end
end

%% Figure
f = figure;
subplot(121);
nseg = [1 2 3 4 5 6 7 10 15 20 30 40 50 52];
plot(nseg, velocity.lotemp.nseg), hold on;
axis([0 55 0 5]);
xticks(0:10:55), xticklabels({'0','10','20','30','40','50',''}), yticks(0:1:5), yticklabels({'0','','','','','5'}), xlabel('# segments'), ylabel('Conduction velocity (m/s)');

subplot(122);
nseg = [1 2 3 4 5 6 7 10 15 20 30 40 50 52];
plot(nseg, velocity.hitemp.nseg), hold on;
axis([0 55 0 5]);
xticks(0:10:55), xticklabels({'0','10','20','30','40','50',''}), yticks(0:1:5), yticklabels({'0','','','','','5'}), xlabel('# segments'), ylabel('Conduction velocity (m/s)');

% Printing as JPEG.
set(gcf,...
    'PaperUnits','inches',...
    'PaperSize',[7 7],...
    'PaperPosition',[0.5 0.5 6 6],...
    'Renderer','Painters');
print(gcf,[saveDirectory '/Cullen2018_R0_Convergence.jpg'],'-djpeg','-r300');
print(gcf,[saveDirectory '/Cullen2018_R0_Convergence.pdf'],'-dpdf');
