clear all;
close all;
% Flag protection through Voronoi Coverage 
% Vincent Maffet 
% 11/26/2018

%
% Variables
%

% Video toggle
videoFLag = 0

% Number of iterations
iterations = 5000;

% Number of robots
N = 12;

% Delta-dist
Delta = 1.5;

% Velocity gain
command_gain = 1;

% Distance for gabriel triangulation
gabriel_dist = 0.3;

% Space element for Voronoi partition
ds = 0.02;

% Voronoi colors
vcolors = rand(N,3);

% Robotarium object
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);

% Area boundaries
bound = r.boundaries;

% Utilities to avoid collisions and change control mapping
si_barrier_cert = create_si_barrier_certificate('SafetyRadius', 1.5*r.robot_diameter);
si_to_uni_dyn = create_si_to_uni_mapping2('LinearVelocityGain', 0.5, ... 
                                          'AngularVelocityLimit', 0.75*r.max_angular_velocity);

% State (consensus > gabriel triangulation > voronoi partition)
state = 'consensus';

% Plot utilities
handlesEdges = cell(N, N);
for i = 1:N
    for j = 1:N
        handlesEdges{i,j} = plot(500, 500, 'black', 'LineWidth', 3);
    end
end

handlesArea = cell((bound(4)-bound(3))/ds, (bound(2)-bound(1))/ds);
i = 1;
for sx = linspace(bound(1), bound(2)-ds, (bound(2)-bound(1))/ds)
    j = 1;
    for sy = linspace(bound(3), bound(4)-ds, (bound(4)-bound(3))/ds)
        handlesArea{j,i} = patch([sx, sx+ds, sx+ds, sx], [sy, sy, sy+ds, sy+ds], 'white', 'LineStyle', 'none');
        handlesArea{j,i}.FaceAlpha = 0 ;
        j = j + 1;
    end
    i = i + 1;
end

handleText = text(0, bound(4)+0.1, state, 'HorizontalAlignment','center','FontSize',14); 

% Flag
plot(0, 0, '.b', 'MarkerSize', 20)

% Initialize video
if videoFLag 
    vid = VideoWriter('flag_protection.avi');
    vid.Quality = 100;
    vid.FrameRate = 72;
    open(vid);
    writeVideo(vid, getframe(gcf));
end

%
% Main loop
%
last = 0;
U = rand(2,N);
for t = 0:iterations
    % Reset edges plot
    for i = 1:N
        for j = 1:N
            handlesEdges{i,j}.XData = 500;
            handlesEdges{i,j}.YData = 500;
        end
    end
    
    X = r.get_poses();
    if max(U(1,:)) < 0.015 && t > last + 50
        last = t
        switch state
            case 'consensus'
                state = 'gabriel';
            case 'gabriel'
                for i = 1:size(handlesArea,2)
                    for j = 1:size(handlesArea,1)
                        handlesArea{j,i}.FaceAlpha = 0.3 ;
                    end
                end
                state = 'voronoi';
            case 'voronoi'
                for i = 1:size(handlesArea,2)
                    for j = 1:size(handlesArea,1)
                        handlesArea{j,i}.FaceAlpha = 0 ;
                    end
                end
                state = 'gabriel';
        end
        handleText.String = state;
    end
    
    U = zeros(2,N);
    
    switch state
        case 'consensus'
            % Agents loop
            for i = 1:N
                U(:, i) = U(:, i) - command_gain*X(1:2, i);
                % Neighbors loop
                for j = delta_disk_neighbors(X, i, Delta)
                    w = 1/(Delta - norm(X(1:2, j) - X(1:2, i)))^2;
                    U(:, i) = U(:, i) + ...
                        command_gain*w*(X(1:2, j) - X(1:2, i));
                end 
            end
            
        case 'gabriel'
            % Agents loop
            for i = 1:N
                
                % Neighbors loop
                for j = delta_disk_neighbors(X, i, Delta)

                    gabriel = true;
                    C = (X(1:2,i)+X(1:2,j))/2;
                    R = norm(X(1:2,i)-C);
                    
                    for k = delta_disk_neighbors(X, i, Delta)
                        if k ~= j && norm(X(1:2,k)-C) < R
                            gabriel = false;
                            break
                        end
                    end

                    if gabriel
                        w = norm(X(1:2, j)-X(1:2, i))^2 - gabriel_dist^2;
                        U(:, i) = U(:, i) + ...
                            command_gain*10*w*(X(1:2, j) - X(1:2, i));
                        handlesEdges{i, j}.XData = [X(1,i); X(1,j)]; 
                        handlesEdges{i, j}.YData = [X(2,i); X(2,j)];
                    end
                end 
            end
            
        case 'voronoi'
            % Space loop
            for sx = linspace(bound(1)+ds/2, bound(2)-ds/2, (bound(2)-bound(1))/ds)
                for sy = linspace(bound(3)+ds/2, bound(4)-ds/2, (bound(4)-bound(3))/ds)
                    sX = [sx ; sy];
                    
                    best = 1;
                    for i = 2:N
                        if norm(sX - X(1:2,i)) < norm(sX - X(1:2,best))
                            best = i;
                        end
                    end
                    
                    U(:,best) = U(:,best) +  ...
                        command_gain*(sX - X(1:2,best))*exp(-norm(sX)^2/0.5);
%                     U(:,best) = U(:,best) +  ...
%                         command_gain*(sX - X(1:2,best))/(2*norm(sX));
                    
                    a = int32(((sx + ds/2) - bound(1))/ds);
                    b = int32(((sy + ds/2) - bound(3))/ds);
                    handlesArea{b,a}.FaceColor = vcolors(best,:);
                end
            end
    end
    
    % Limit actuator saturation
    for i = 1:N
        if norm(U(:,i)) > r.max_linear_velocity
            factor = r.max_linear_velocity / norm(U(:,i));
            U(:,i) = factor * U(:,i);
        end   
    end
    % Avoid boundaries
    dt = r.time_step;
    for i = 1:N
        if ((X(1,i)+dt*U(1,i)) < bound(1) || (X(1,i)+dt*U(1,i)) > bound(2))
            U(1,i) = 0;
        end
        if ((X(2,i)+dt*U(2,i)) < bound(3) || (X(2,i)+dt*U(2,i)) > bound(4))
            U(2,i) = 0;
        end
    end
    % Avoid collisions
    U = si_barrier_cert(U, X);
    % Transform the single-integrator dynamics to unicycle dynamics
    U = si_to_uni_dyn(U, X);
    
    % Apply command
    r.set_velocities(1:N, U);
    r.step();
    
    if videoFLag && mod(t,5)                               % Record a video frame every 10 iterations
            writeVideo(vid, getframe(gcf)); 
    end
    
end

if videoFLag; close(vid); end

r.debug();