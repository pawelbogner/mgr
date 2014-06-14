results=zeros(10, 7);
counter=1;
for longitudinal=[1 15]
    for lateral=[1 15]
        for time_horizon=[20 10]
            basename=['final_' num2str(longitudinal) '_' num2str(lateral) '_' num2str(time_horizon) ]
            name=[basename '.mat' ];
            if exist(name, 'file')
                load(name);
                figure(1);
                set(gca,'FontSize',17)
                plot(x_ksi(:,1), x_ksi(:,2), 'k');
                set(gca,'FontSize',17)
                xlabel('x');
                ylabel('y');
                print(1, [basename '_path'], '-depsc')
                print(1, [basename '_path'], '-dpng')
                figure(2);
                set(gca,'FontSize',17)
                s=calc_slips(x_ksi(:,1:10));
                subplot(2, 2, 1);
                set(gca,'FontSize',17)
                plot(t, s(:,1), 'k');
                xlabel('t');
                ylabel('s_{14}');
                subplot(2, 2, 2);
                set(gca,'FontSize',17)
                plot(t, s(:,2), 'k');
                xlabel('t');
                ylabel('s_{23}');
                subplot(2, 2, 3);
                set(gca,'FontSize',17)
                plot(t, s(:,3), 'k');
                xlabel('t');
                ylabel('s_{12}');
                subplot(2, 2, 4);
                set(gca,'FontSize',17)
                plot(t, s(:,4), 'k');
                xlabel('t');
                ylabel('s_{34}');
                print(2, [basename '_slips'], '-depsc')
                print(2, [basename '_slips'], '-dpng')
                figure(3);
                u=calc_u(t, lambda, time_horizon);
                plot(t, u(:,1), 'k-', t, u(:,2), 'k--');
                set(gca,'FontSize',17)
                ylabel('u');
                xlabel('t');
                legend('u_1', 'u_2', 'Location', 'best');
                print(3, [basename '_u'], '-depsc')
                print(3, [basename '_u'], '-dpng')
                results(counter,:)=[longitudinal lateral time_horizon trapz(u.^2) max(abs(u))]
                counter=counter+1;
            end
        end
    end
end


basename='manip_task'
name=[basename '.mat' ];
if exist(name, 'file')
    load(name);
    figure(1);
    plot(x_ksi(:,1), x_ksi(:,2), 'k');
                set(gca,'FontSize',17)
    xlabel('x','FontSize',17);
    ylabel('y','FontSize',17);
    print(1, [basename '_path_manip'], '-depsc')
    print(1, [basename '_path_manip'], '-dpng')
    figure(2);
                set(gca,'FontSize',17)
    s=calc_slips(x_ksi(:,1:10));
    subplot(2, 2, 1);
                set(gca,'FontSize',17)
    plot(t, s(:,1), 'k');
    xlabel('t','FontSize',17);
    ylabel('s14','FontSize',17);
    subplot(2, 2, 2);
                set(gca,'FontSize',17)
    plot(t, s(:,2), 'k');
    xlabel('t','FontSize',17);
    ylabel('s23','FontSize',17);
    subplot(2, 2, 3);
                set(gca,'FontSize',17)
    plot(t, s(:,3), 'k');
    xlabel('t','FontSize',17);
    ylabel('s12','FontSize',17);
    subplot(2, 2, 4);
                set(gca,'FontSize',17)
    plot(t, s(:,4), 'k');
    xlabel('t','FontSize',17);
    ylabel('s34','FontSize',17);
    print(2, [basename '_slips_manip'], '-depsc')
    print(2, [basename '_slips_manip'], '-dpng')
    figure(3);
    u=calc_u(t, lambda, 20);
                plot(t, u(:,1), 'k-', t, u(:,2), 'k--');
                set(gca,'FontSize',17)
    ylabel('u','FontSize',17);
    xlabel('t','FontSize',17);
    legend('u_1', 'u_2', 'Location', 'best');
    print(3, [basename '_u_manip'], '-depsc')
    print(3, [basename '_u_manip'], '-dpng')
    results(counter,:)=[15 1 20 trapz(u.^2) max(abs(u))]
    counter=counter+1;
end

basename='manip_pltf_task'
name=[basename '.mat' ];
if exist(name, 'file')
    load(name);
    figure(1);
    plot(x_ksi(:,1), x_ksi(:,2), 'k');
                set(gca,'FontSize',17)
    xlabel('x','FontSize',17);
    ylabel('y','FontSize',17);
    print(1, [basename '_path_manip_pltf'], '-depsc')
    print(1, [basename '_path_manip_pltf'], '-dpng')
    figure(2);
    s=calc_slips(x_ksi(:,1:10));
    subplot(2, 2, 1);
                set(gca,'FontSize',17)
    plot(t, s(:,1), 'k');
    xlabel('t','FontSize',17);
    ylabel('s14','FontSize',17);
    subplot(2, 2, 2);
                set(gca,'FontSize',17)
    plot(t, s(:,2), 'k');
    xlabel('t','FontSize',17);
    ylabel('s23','FontSize',17);
    subplot(2, 2, 3);
                set(gca,'FontSize',17)
    plot(t, s(:,3), 'k');
    xlabel('t','FontSize',17);
    ylabel('s12','FontSize',17);
    subplot(2, 2, 4);
                set(gca,'FontSize',17)
    plot(t, s(:,4), 'k');
    xlabel('t','FontSize',17);
    ylabel('s34','FontSize',17);
    print(2, [basename '_slips_manip_pltf'], '-depsc')
    print(2, [basename '_slips_manip_pltf'], '-dpng')
    figure(3);
    u=calc_u(t, lambda, 20);
                plot(t, u(:,1), 'k-', t, u(:,2), 'k--');
                set(gca,'FontSize',17)
    ylabel('u','FontSize',17);
    xlabel('t','FontSize',17);
    legend('u_1', 'u_2', 'Location', 'best');
    print(3, [basename '_u_manip_pltf'], '-depsc')
    print(3, [basename '_u_manip_pltf'], '-dpng')
    results(counter,:)=[15 1 20 trapz(u.^2) max(abs(u))]
    counter=counter+1;
end

save('energy_ampl.mat', 'results')
