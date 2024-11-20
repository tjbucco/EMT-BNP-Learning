function plotMovingTargets(data_in, hyper)

% plot generated data

Color = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
figure
%title('Target tracking')
axis equal
hold on
%set(0, 'DefaultLineLineWidth', 3);
%set(0, 'DefaultLineMarkerSize', 6);
for t = 1:data_in.maxTactive
    for n = 1:data_in.maxNactive
        if t == data_in.Ts(n) 
            %plot(data_in.X(1, 1, n), data_in.X(2, 1, n), '+-', 'Color',Color{mod(n - 1, 5) + 1}, "LineWidth",1,"MarkerSize",6)
            plot(data_in.X(1, t, n), data_in.X(2, t, n), "diamond", 'Color','k', "LineWidth",1,"MarkerSize",6)
            %text(data_in.X(1, t, n) + 50*(-1)^(n), data_in.X(2, t, n) + 50*(-1)^(n), sprintf("k = %i", round(t)));
        elseif any(t == data_in.Ts(n):data_in.Te(n))
            plot(data_in.X(1, (t - 1):t, n), data_in.X(2, (t - 1):t, n), '+-', 'Color',Color{mod(data_in.C(n), 5) + 1}, "LineWidth",0.5,"MarkerSize",6)
        end
        egtnr = calcelipse(data_in.X(:, t, n), data_in.Dstar(:, :, data_in.C(n)));
        plot(egtnr(1,:), egtnr(2,:), 'Color','k',"LineWidth",1);
        if t == data_in.maxTactive
            plot(data_in.Z(1, :, t, n), data_in.Z(2, :, t, n), '.', 'Color',Color{mod(data_in.C(n), 5) + 1},"MarkerSize",10);
        end
    end
    
    if size(data_in.Z, 4) > data_in.maxNactive
        plot(data_in.Z(1, :, t, data_in.maxNactive+1), data_in.Z(2, :, t, data_in.maxNactive+1), '.', 'Color','k',"MarkerSize",10);
    end
end
for t = 1:data_in.maxTactive
    for n = 1:data_in.maxNactive
        if t == data_in.Ts(n)
            %plot(data_in.X(1, 1, n), data_in.X(2, 1, n), '+-', 'Color',Color{mod(n - 1, 5) + 1}, "LineWidth",1,"MarkerSize",6)
            plot(data_in.X(1, t, n), data_in.X(2, t, n), "diamond", 'Color','k', "LineWidth",2,"MarkerSize",7)
            %plot(data_in.X(1, t, n), data_in.X(2, t, n), "diamond", 'Color','k', "MarkerFaceColor", "#FF0000", "MarkerEdgeColor","k", "LineWidth",2,"MarkerSize",7)
            %text(data_in.X(1, t, n) + 50*(-1)^(n), data_i
        end
    end
end

set(gca, 'FontName', 'Times')
set(gca, 'fontsize', 12)
set(0, 'DefaultTextFontName', 'Times')
set(0, 'DefaultTextFontSize', 12)
% text(data_in.X(1, data_in.Ts(1), 1) + 1, data_in.X(2, data_in.Ts(1), 1) + 110, sprintf("$k = %i$", round(data_in.Ts(1))), "Interpreter","latex");
% text(data_in.X(1, data_in.Ts(2), 2) -30, data_in.X(2, data_in.Ts(2), 2) + 100, sprintf("$k = %i$", round(data_in.Ts(2))), "Interpreter","latex");
% text(data_in.X(1, data_in.Ts(3), 3) + 1, data_in.X(2, data_in.Ts(3), 3) - 110, sprintf("$k = %i$", round(data_in.Ts(3))), "Interpreter","latex");
% text(data_in.X(1, data_in.Ts(4), 4) - 50, data_in.X(2, data_in.Ts(4), 4) + 100, sprintf("$k = %i$", round(data_in.Ts(4))), "Interpreter","latex");
% text(data_in.X(1, data_in.Ts(5), 5) - 50, data_in.X(2, data_in.Ts(5), 5) - 110, sprintf("$k = %i$", round(data_in.Ts(5))), "Interpreter","latex");
% text(data_in.X(1, data_in.Ts(6), 6) - 150, data_in.X(2, data_in.Ts(6), 6) + 110, sprintf("$k = %i$", round(data_in.Ts(6))), "Interpreter","latex");


%axis([1.1*data_in.min_xy(1), 1.1*data_in.max_xy(1), 1.1*data_in.min_xy(2), 1.1*data_in.max_xy(2)]);
xlim([-hyper.radius-50, hyper.radius+50]); ylim([-hyper.radius-50, hyper.radius+50])
xticks(linspace(-hyper.radius,hyper.radius,6)); yticks(linspace(-hyper.radius,hyper.radius,6));
xlabel("Dimension 1 [m]"); ylabel("Dimension 2 [m]");
%legend({'x.^1';'x.^2';'x.^3';'x.^4';'x.^5';'x.^6';'x.^7'})

%pause;

end
