close all
clear all

models = {'dokos_model_1996'};
blocks = 0:10:100;
chaste_test_output_path = '/export/testoutput';

for m=1:length(models)
    
    apds = zeros(1:length(blocks),1);
    
    for i=1:length(blocks)
        
        d = importdata([chaste_test_output_path '/SanWithFunnyCurrentBlock/' models{m} '_block_' num2str(blocks(i))]);
        
        t = d(:,1);
        
        % Get the gradients of voltage
        v_dot = diff(d(:,2));
        v_double_dot = diff(v_dot);
        
        % Get somewhere where the gradient is negative
        idx = intersect(find(v_dot < 0), find(t>500));
        
        % At the end of the first continuous run of negative dt is a minimum of the voltage.
        for j=1:length(idx)
            expected = idx(1) + j - 1;
            if expected ~= idx(j)
                minimum_at_t = idx(j-1);
                break
            end
        end
        
        % After this, find the first time voltage goes above APD50
        % threshold
        threshold = min(d(:,2)) + 0.5*(max(d(:,2)) - min(d(:,2)));
        idx_above = minimum_at_t + find(d(minimum_at_t:end,2)>threshold, 1, 'first');
        time_above = t(idx_above);
        time_below = t(idx_above + find(d(idx_above:end,2)<threshold, 1, 'first'));
        
        apds(i) = time_below - time_above;
        
        % I just used this to get the minimum detection working nicely..!
        % figure(2)
        % plot(t(idx),zeros(length(idx),1),'r.')
        % hold on
        % plot(t(1:end-1), v_dot)
        
        figure(m)
        subplot(1,2,1)
        plot(t(minimum_at_t:end) - t(minimum_at_t), d(minimum_at_t:end,2), '-')
        hold all
        
    end
    title(models{m})
    xlim([0 450])
    ylim([-75 30])
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    
    subplot(1,2,2)
    plot(blocks,100.0*(apds./apds(1))-100.0,'b.-')
    xlim([0 100])
    xlabel('Funny Current Block (%)')
    ylabel('Prolongation of APD_{50} (%)')
end
