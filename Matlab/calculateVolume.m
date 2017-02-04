% Christina Aiello, 11/8/2015
% This function will calculate the volume of data 
% (meant for gyroscope data).

function volume = calculateVolume(gyroX, gyroY, gyroZ, draw);
    % Reference:
    % https://www.mathworks.com/matlabcentral/answers/253258-plot-a-circle-around-xyz-plotted-points-and-then-get-its-volume

    if (nargin == 3)
        draw = false;
    end
    xyz = [gyroX(:) gyroY(:) gyroZ(:)];                     % Consolidated: [xVector(:) yVector(:) zVector(:)]
    C = mean(xyz);                                          % Calculate Centre
    R = sum(bsxfun(@minus, xyz, C).^2, 2);                  % Calculate Radii
    [Rmax, idx] = max(R);                                   % Calculate Maximum Radius
    volume = 4*pi*Rmax^3/3;                                 % Volume Of Sphere
    if(draw)
        [Xs,Ys,Zs] = sphere;                                    % Use sphere Function
        Rsph = sqrt(Rmax);                                      % Calculate Radius Of Sphere
        figure(1);
        plot3(xyz(:,1), xyz(:,2), xyz(:,3), '.');                % Plot Data
        hold on;
        plot3(C(1), C(2), C(3), 'r*');                           % Plot Centre
        plot3(xyz(idx,1), xyz(idx,2), xyz(idx,3), 'bp');         % Plot Maximum Radius
        Sph = mesh(Rsph*Xs+C(1), Rsph*Ys+C(2), Rsph*Zs+C(3));   % Calculate Enveloping Sphere
        hold off;
        grid on;
        axis equal;
        set(Sph, 'FaceAlpha',0);                               % Set Sphere Transparency
        set(Sph, 'EdgeAlpha',0.4);                               % Set Sphere Transparency
        legend('Data', 'Centre', 'Max Radius', 'Location','EastOutside');
    end
    % Step 1: Concatenating the term 'Sway Area,' the file name for this
    % trial, and the trial number:
    % Lastly, clear the figure:
    
end