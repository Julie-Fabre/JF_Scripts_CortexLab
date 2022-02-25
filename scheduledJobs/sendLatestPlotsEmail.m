
%For this to work in Gmail you must allow less secure connections in gmail.

%% generate plots 
close all;
errors = '';
thisDate = datestr(datetime('yesterday', 'Format', 'yyyy-MM-dd'), 'yyyy-mm-dd');
% run KS
try
    JF_preprocess_NPX2('JF067', thisDate, 'neuropixPhase3B1_kilosortChanMap.mat', 1, 1, [],1);
catch
    errors = [errors, 'error KSing '];
end
% get behavior + save plots
try
    close all;
    if thisDate > datetime(2022,02,26) %Laura starts training 
        bhvData = noGoWorld_behavior({'JF067', 'JF070', 'JF071', 'JF072', 'JF073'}); %behavior plots
        saveas(figure(6), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorLines.png'])%save behavior plots 
        saveas(figure(7), ['~/Documents/cronJobPlots/', thisDate, '_JF070_behaviorLines.png'])%save behavior plots 
        saveas(figure(8), ['~/Documents/cronJobPlots/', thisDate, '_JF071_behaviorLines.png'])%save behavior plots 
        saveas(figure(9), ['~/Documents/cronJobPlots/', thisDate, '_JF072_behaviorLines.png'])%save behavior plots 
        saveas(figure(10), ['~/Documents/cronJobPlots/', thisDate, '_JF073_behaviorLines.png'])%save behavior plots 
        screensize = get( groot, 'Screensize' );
        set(figure(1), 'Position',  screensize)
        set(figure(2), 'Position',  screensize)
        set(figure(3), 'Position',  screensize)
        set(figure(4), 'Position',  screensize)
        set(figure(5), 'Position',  screensize)
        saveas(figure(1), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorDetails.png'])%save behavior plots
        saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF070_behaviorDetails.png'])%save behavior plots
        saveas(figure(3), ['~/Documents/cronJobPlots/', thisDate, '_JF071_behaviorDetails.png'])%save behavior plots
        saveas(figure(4), ['~/Documents/cronJobPlots/', thisDate, '_JF072_behaviorDetails.png'])%save behavior plots
        saveas(figure(5), ['~/Documents/cronJobPlots/', thisDate, '_JF073_behaviorDetails.png'])%save behavior plots
        close all;
    else
        bhvData = noGoWorld_behavior({'JF067'}); %behavior plots
        saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorLines.png'])%save behavior plots 
        screensize = get( groot, 'Screensize' );
        set(figure(1), 'Position',  screensize)
        saveas(figure(1), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorDetails.png'])%save behavior plots 
        close all;
    end
    
catch
    errors = [errors, 'error behavior '];
end
try
    chronicOverLearning;% get chronic over learning plots 
    saveas(figure(3), ['~/Documents/cronJobPlots/', thisDate, '_JF067_driftMap.png'])%save behavior plots 
    saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF067_unitNumbers.png'])%save behavior plots 
    screensize = get( groot, 'Screensize' );
    set(figure(4), 'Position',  screensize)
    saveas(figure(4), ['~/Documents/cronJobPlots/', thisDate, '_JF067_PSTH.png'])%save behavior plots 
    close all;
catch
    errors = [errors, 'error PSTH '];
end

%% load email credentials
credentials = readtable('~/credentials.txt');

%% setup SMTP
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
setpref('Internet', 'E_mail', credentials.email{:});
setpref('Internet', 'SMTP_Username', credentials.email{:}(1:end-10));
setpref('Internet', 'SMTP_Password', credentials.password{:});

%% Gmail server
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth', 'true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port', '465');

%% fetch last plots 
plots = dir(['~/Documents/cronJobPlots/' thisDate '*JF067*.png']);
theseAttachedPlots = cell(length(plots),1);
for iPlot = 1:length(plots)
    theseAttachedPlots{iPlot} = fullfile(plots(iPlot).folder, plots(iPlot).name);
end

if thisDate > datetime(2022,02,26)
    plotsAcuteMice = dir(['~/Documents/cronJobPlots/' thisDate '*JF07*.png']);
    theseAttachedAcutePlots = cell(length(plotsAcuteMice),1);
    for iPlot = 1:length(theseAttachedAcutePlots)
        theseAttachedAcutePlots{iPlot} = fullfile(plotsAcuteMice(iPlot).folder, plotsAcuteMice(iPlot).name);
    end
end

%% send plots by email 
message = ['Hey there bb, here are your behavior and PSTH plots for ', thisDate '.'];
sendmail(credentials.mailto{:}, [thisDate ' JF067 plots '], [message, ' ', errors], theseAttachedPlots)
sendmail(credentials.mailto{:}, [thisDate ' acute mice training plots '], [message, ' ', errors], theseAttachedAcutePlots)

%%
exit; %close matlab


























