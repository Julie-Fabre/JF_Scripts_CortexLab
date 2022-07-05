
%For this to work in Gmail you must allow less secure connections in gmail.

%% generate plots 
close all;
errors = '';
thisDate = datestr(datetime('today', 'Format', 'yyyy-MM-dd'), 'yyyy-mm-dd');
% run KS
try
    JF_preprocess_NPX2('JF067', thisDate, 'neuropixPhase3B1_kilosortChanMap.mat', 1, 1, [],1);
    JF_preprocess_NPX2('JF071', thisDate, 'neuropixPhase3B1_kilosortChanMap.mat', 1, 1, [],1);
    JF_preprocess_NPX2('JF071', thisDate, 'neuropixPhase3B1_kilosortChanMap.mat', 1, 2, [],1);
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
        pause(50)
        saveas(figure(1), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorDetails.png'])%save behavior plots
        saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF070_behaviorDetails.png'])%save behavior plots
        saveas(figure(3), ['~/Documents/cronJobPlots/', thisDate, '_JF071_behaviorDetails.png'])%save behavior plots
        saveas(figure(4), ['~/Documents/cronJobPlots/', thisDate, '_JF072_behaviorDetails.png'])%save behavior plots
        saveas(figure(5), ['~/Documents/cronJobPlots/', thisDate, '_JF073_behaviorDetails.png'])%save behavior plots
        close all;
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF067_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/83c1549c-8ad4-42d5-adf7-0ca9280a5070');
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF070_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/04d9a0c3-01cc-4ca4-86c2-e48da5284333');
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF071_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/12119bed-6e0e-4b8e-94d2-ed8186430261');
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF072_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/4e6b318e-990a-48cb-8e29-fc7f3cffde64');
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF063_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/fd301c93-cd90-46af-8629-80ce6340c5c7');
    
    else
        bhvData = noGoWorld_behavior({'JF067'}); %behavior plots
        saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorLines.png'])%save behavior plots 
        screensize = get( groot, 'Screensize' );
        set(figure(1), 'Position',  screensize)
        pause(10)
        saveas(figure(1), ['~/Documents/cronJobPlots/', thisDate, '_JF067_behaviorDetails.png'])%save behavior plots 
        close all;
        websave(['~/Documents/cronJobPlots/', thisDate, '_JF067_weight.html'],'https://alyx.cortexlab.net/admin-actions/water-history/83c1549c-8ad4-42d5-adf7-0ca9280a5070');
  
    end
    
catch
    errors = [errors, 'error behavior '];
end
% get PSTH
try
    chronicOverLearning;% get chronic over learning plots 
    saveas(figure(3), ['~/Documents/cronJobPlots/', thisDate, '_JF067_driftMap.png'])%save behavior plots 
    saveas(figure(2), ['~/Documents/cronJobPlots/', thisDate, '_JF067_unitNumbers.png'])%save behavior plots 
    screensize = get( groot, 'Screensize' );
    set(figure(4), 'Position',  screensize)
    pause(30)
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
plots = dir(['~/Documents/cronJobPlots/' thisDate '*JF067*']);
theseAttachedPlots = cell(length(plots),1);
for iPlot = 1:length(plots)
    theseAttachedPlots{iPlot} = fullfile(plots(iPlot).folder, plots(iPlot).name);
end

if thisDate > datetime(2022,02,26)
    plotsAcuteMice = dir(['~/Documents/cronJobPlots/' thisDate '*JF07*']);
    theseAttachedAcutePlots = cell(length(plotsAcuteMice),1);
    for iPlot = 1:length(theseAttachedAcutePlots)
        theseAttachedAcutePlots{iPlot} = fullfile(plotsAcuteMice(iPlot).folder, plotsAcuteMice(iPlot).name);
    end
end

%% get a joke 
try 
   
    contents = urlread(' https://icanhazdadjoke.com/');
    searchString = '<meta property="og:description" content="';
    location = strfind(contents, searchString);
    locationEnd = strfind(contents(location:end), '" />');
    joke = contents(location+length(searchString):location+locationEnd(1)-2);
    
    contents2 = urlread(' https://icanhazdadjoke.com/');
    location2 = strfind(contents2, searchString);
    locationEnd2 = strfind(contents2(location2:end), '" />');
    joke2 = contents2(location2+length(searchString):location2+locationEnd2(1)-2);
catch
    joke = '';
    joke2 = '';
end
%% send plots by email 
message = ['Hey Julie, ' newline joke newline 'Here are your behavior and PSTH plots for ', thisDate '.' newline errors newline 'Cheers,' newline 'Automated Julie'];
message2 = ['Hey Julie, ' newline joke2 newline 'Here are your behavior plots for ', thisDate '.' newline errors newline 'Cheers,' newline 'Automated Julie'];
sendmail(credentials.mailto{:}, [thisDate ' JF067 plots '], message, theseAttachedPlots)
sendmail(credentials.mailto{:}, [thisDate ' acute mice training plots '], message2, theseAttachedAcutePlots)

%%
exit; %close matlab

























