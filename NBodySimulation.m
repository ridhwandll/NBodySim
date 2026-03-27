classdef NBodySimulation < matlab.apps.AppBase
   
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        MainGrid        matlab.ui.container.GridLayout
        ControlPanel    matlab.ui.container.Panel
        CtrlGrid        matlab.ui.container.GridLayout
        TitleLabel      matlab.ui.control.Label
        NumBodiesSlider matlab.ui.control.Slider
        NumBodiesValue  matlab.ui.control.Label
        PresetDropDown  matlab.ui.control.DropDown
        GravSlider      matlab.ui.control.Slider
        SoftenSlider    matlab.ui.control.Slider
        TimeStepSlider  matlab.ui.control.Slider
        TrailSlider     matlab.ui.control.Slider
        StartButton     matlab.ui.control.Button
        StopButton      matlab.ui.control.Button
        ResetButton     matlab.ui.control.Button
        StatsPanel      matlab.ui.container.Panel
        StatsGrid       matlab.ui.container.GridLayout
        KineticValue    matlab.ui.control.Label
        SimTimeValue    matlab.ui.control.Label
        StepValue       matlab.ui.control.Label
        SimAxes         matlab.ui.control.UIAxes
    end

    properties (Access = private)
        N, pos, vel, mass, colors, radii
        trailX, trailY, trailLen
        bodyPlots, trailLines
        G, dt, softening, simTime, stepCount
        simTimer, isRunning = false, isInit = false
    end

    methods (Access = private)
        function initSimulation(app)
            app.stopSim(); % Ensure timer is stopped before re-init
            app.N         = round(app.NumBodiesSlider.Value);
            app.G         = app.GravSlider.Value;
            app.dt        = app.TimeStepSlider.Value;
            app.softening = app.SoftenSlider.Value;
            app.trailLen  = round(app.TrailSlider.Value);
            app.simTime   = 0;
            app.stepCount = 0;
            app.setupPreset(app.PresetDropDown.Value);
            app.initGraphics();
            app.isInit = true;
            app.updateStats();
        end

        function setupPreset(app, preset)
            N = app.N; G = app.G;
            app.mass = ones(N,1); 
            switch preset
                case 'Random Cluster'
                    R = 25;
                    app.pos = (rand(N,2)-0.5) * 2 * R;
                    vScale = sqrt(G * N / R) * 0.4;
                    app.vel = (rand(N,2)-0.5) * vScale;
                case 'Solar System'
                    N = min(N,12); app.N = N;
                    app.mass = ones(N,1) * 0.2;
                    app.mass(1) = 150; 
                    app.pos = zeros(N,2); app.vel = zeros(N,2);
                    for i = 2:N
                        r = 12 + (i-2)*7;
                        theta = rand() * 2 * pi;
                        app.pos(i,:) = [r*cos(theta), r*sin(theta)];
                        v = sqrt(G * app.mass(1) / r);
                        app.vel(i,:) = [-sin(theta), cos(theta)] * v;
                    end
                case 'Binary System'
                    m_star = 80; sep = 12;
                    app.mass(1:2) = m_star;
                    app.pos(1,:) = [-sep, 0]; app.pos(2,:) = [sep, 0];
                    v_orbit = sqrt(G * m_star / (4 * sep));
                    app.vel(1,:) = [0, v_orbit]; app.vel(2,:) = [0, -v_orbit];
                    for i = 3:N
                        r = 35 + rand()*15; th = rand()*2*pi;
                        app.pos(i,:) = [r*cos(th), r*sin(th)];
                        app.vel(i,:) = [-sin(th), cos(th)] * sqrt(G*2*m_star/r);
                    end
            end
            com_v = sum(app.mass .* app.vel) / sum(app.mass);
            app.vel = app.vel - com_v;
            app.colors = hsv(app.N);
            app.radii  = (5 + 12 * (app.mass/max(app.mass)).^(1/3)).^2;
            app.trailX = cell(N,1); app.trailY = cell(N,1);
            for i=1:N, app.trailX{i}=app.pos(i,1); app.trailY{i}=app.pos(i,2); end
        end

        function initGraphics(app)
            cla(app.SimAxes); hold(app.SimAxes, 'on');
            sPos = (rand(200,2)-0.5)*1200;
            scatter(app.SimAxes, sPos(:,1), sPos(:,2), 1.5, [0.3 0.3 0.5], 'filled', 'HandleVisibility', 'off');
            app.trailLines = gobjects(app.N,1);
            for i = 1:app.N
                app.trailLines(i) = line(app.SimAxes, app.pos(i,1), app.pos(i,2), ...
                    'Color', [app.colors(i,:), 0.4], 'LineWidth', 1.2);
            end
            app.bodyPlots = scatter(app.SimAxes, app.pos(:,1), app.pos(:,2), ...
                app.radii, app.colors, 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.5);
            hold(app.SimAxes, 'off');
        end

        function stepSimulation(app)
            m = app.mass; G = app.G; dt = app.dt; e2 = app.softening^2;
            a1 = nbodyAcc(app.pos, m, G, e2);
            v_mid = app.vel + a1 * (dt/2);
            app.pos = app.pos + v_mid * dt;
            a2 = nbodyAcc(app.pos, m, G, e2);
            app.vel = v_mid + a2 * (dt/2);
            app.simTime = app.simTime + dt;
            app.stepCount = app.stepCount + 1;
        end

        function updateGraphics(app)
            if ~isvalid(app.UIFigure), return; end
            for i = 1:app.N
                app.trailX{i}(end+1) = app.pos(i,1);
                app.trailY{i}(end+1) = app.pos(i,2);
                if numel(app.trailX{i}) > app.trailLen
                    app.trailX{i}(1) = []; app.trailY{i}(1) = [];
                end
                set(app.trailLines(i), 'XData', app.trailX{i}, 'YData', app.trailY{i});
            end
            set(app.bodyPlots, 'XData', app.pos(:,1), 'YData', app.pos(:,2));
            limit = max(max(abs(app.pos))) * 1.3 + 10;
            axis(app.SimAxes, [-limit limit -limit limit]);
            drawnow limitrate;
        end

        function updateStats(app)
            ke = 0.5 * sum(app.mass .* sum(app.vel.^2, 2));
            app.KineticValue.Text = sprintf('KE: %.2f', ke);
            app.SimTimeValue.Text = sprintf('Time: %.2f', app.simTime);
            app.StepValue.Text    = sprintf('Steps: %d', app.stepCount);
        end

        function timerFcn(app, ~, ~)
            if ~app.isRunning || ~isvalid(app.UIFigure), return; end
            try
                for k = 1:4
                    app.stepSimulation();
                end
                app.updateGraphics();
                if mod(app.stepCount, 40) == 0, app.updateStats(); end
            catch
                app.stopSim();
            end
        end

        function createComponents(app)
            bg_dark = [0.01 0.01 0.01];
            panel_dark = [0.01 0.01 0.11];
            text_color = [0.8 0.8 0.8];
            accent = [0.8 0.8 0.8];

            app.UIFigure = uifigure('Name','N-Body Deep Space','Color', bg_dark);
            app.UIFigure.Position = [100 100 1100 700];
            app.MainGrid = uigridlayout(app.UIFigure, [1 2], 'ColumnWidth', {280, '1x'}, 'BackgroundColor', bg_dark);
            
            app.ControlPanel = uipanel(app.MainGrid, 'BackgroundColor', panel_dark, 'BorderType', 'none');
            app.CtrlGrid = uigridlayout(app.ControlPanel, [14 1], 'RowHeight', repmat({30}, 1, 14), 'BackgroundColor', panel_dark);
            
            app.TitleLabel = uilabel(app.CtrlGrid, 'Text', 'GRAVITY SIMULATOR', 'FontSize', 16, 'FontWeight', 'bold', 'FontColor', accent);

            %TODO FIX THIS ___________________
            uilabel(app.CtrlGrid, 'Text', 'Gravity Constant', 'FontColor', text_color);
            app.GravSlider = uislider(app.CtrlGrid, 'Limits', [0.1 5], 'Value', 1);
            uilabel(app.CtrlGrid, 'Text', 'Temporal Step (dt)', 'FontColor', text_color);
            app.TimeStepSlider = uislider(app.CtrlGrid, 'Limits', [0.005 0.1], 'Value', 0.02);
            uilabel(app.CtrlGrid, 'Text', 'Softening Radius', 'FontColor', text_color);
            app.SoftenSlider = uislider(app.CtrlGrid, 'Limits', [0.1 5], 'Value', 0.8);
            app.PresetDropDown = uidropdown(app.CtrlGrid, 'Items', {'Random Cluster', 'Solar System', 'Binary System'});
            uilabel(app.CtrlGrid, 'Text', 'Object Count', 'FontColor', text_color);
            app.NumBodiesSlider = uislider(app.CtrlGrid, 'Limits', [2 60], 'Value', 15);
            uilabel(app.CtrlGrid, 'Text', 'Trail Persistance', 'FontColor', text_color);
            app.TrailSlider = uislider(app.CtrlGrid, 'Limits', [5 1000], 'Value', 150);
            %TODO FIX THIS ^^^^^^^^^^^^^^^^^^^

            app.StartButton = uibutton(app.CtrlGrid, 'Text', 'PLAY / RESUME', 'BackgroundColor', [0.1 0.4 0.2], 'FontColor', 'w', 'ButtonPushedFcn', @(~,~) app.startSim());
            app.StopButton  = uibutton(app.CtrlGrid, 'Text', 'PAUSE', 'BackgroundColor', [0.4 0.1 0.1], 'FontColor', 'w', 'ButtonPushedFcn', @(~,~) app.stopSim());
            app.ResetButton = uibutton(app.CtrlGrid, 'Text', 'RESET SYSTEM', 'BackgroundColor', [0.2 0.2 0.4], 'FontColor', 'w', 'ButtonPushedFcn', @(~,~) app.initSimulation());

            app.StatsPanel = uipanel(app.CtrlGrid, 'Title', 'System Metrics', 'BackgroundColor', panel_dark, 'ForegroundColor', accent);
            app.StatsGrid = uigridlayout(app.StatsPanel, [3 1], 'RowHeight', {20, 20, 20}, 'Padding', [5 5 5 5], 'BackgroundColor', panel_dark);
            app.KineticValue = uilabel(app.StatsGrid, 'Text', 'KE: 0', 'FontColor', [0.6 1.0 0.6]);
            app.SimTimeValue = uilabel(app.StatsGrid, 'Text', 'Time: 0', 'FontColor', [1.0 0.8 0.4]);
            app.StepValue = uilabel(app.StatsGrid, 'Text', 'Steps: 0', 'FontColor', [0.4 0.8 1.0]);

            app.SimAxes = uiaxes(app.MainGrid, 'BackgroundColor', bg_dark, 'XColor', [0.2 0.2 0.3], 'YColor', [0.2 0.2 0.3]);
            app.SimAxes.Toolbar.Visible = 'off';
            app.SimAxes.Color = [0 0 0];
            app.SimAxes.XTickLabel = {};
            app.SimAxes.YTickLabel = {};
            grid(app.SimAxes, 'off');
            app.SimAxes.GridAlpha = 0.1;
        end

        function startSim(app)
            if ~app.isInit, app.initSimulation(); end
            app.isRunning = true;
            % The fix: only start if it's currently stopped
            if isvalid(app.simTimer) && strcmp(app.simTimer.Running, 'off')
                start(app.simTimer);
            end
        end

        function stopSim(app)
            app.isRunning = false;
            if isvalid(app.simTimer) && strcmp(app.simTimer.Running, 'on')
                stop(app.simTimer);
            end
        end
    end

    methods (Access = public)
        function app = NBodySimulation
            createComponents(app);
            app.simTimer = timer('ExecutionMode','fixedRate','Period',0.033,'TimerFcn',@app.timerFcn, 'BusyMode', 'drop');
            app.initSimulation();
        end
        function delete(app)
            app.stopSim();
            if isvalid(app.simTimer), delete(app.simTimer); end
        end
    end
end

function acc = nbodyAcc(pos, m, G, e2)
    N = size(pos,1);
    dx = pos(:,1)' - pos(:,1); 
    dy = pos(:,2)' - pos(:,2);
    r2 = dx.^2 + dy.^2 + e2;
    invR3 = r2.^(-1.5);
    invR3(1:N+1:end) = 0; 
    fac = G * (m') .* invR3;
    acc = [sum(fac.*dx, 2), sum(fac.*dy, 2)];
end