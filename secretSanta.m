


%% Email Preferences
%For this to work in Gmail you must allow less secure connections in gmail.
setpref('Internet','SMTP_Server','smtp.gmail.com');
%Your email goes below.
setpref('Internet','E_mail','username@gmail.com');
%Your username.
setpref('Internet','SMTP_Username','username');
%Your password.
setpref('Internet','SMTP_Password','password');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');



%% People being loaded in.
people = {};
test = 0; %If test = 1, will send out a test email to everyone. If test = 0,
% will send actual secret santa email.

%The first column in the matrix corresponds to someones name. The second
%corresponds to their email.
people{1,1} = 'Julie';
people{1,2} = 'fakeemail4@gmail.com';


%Finding the number of people.
numPeople = size(people,1);
%A simple matrix to be used for comparison. The Gift Giver
numHolder = linspace(1, numPeople, numPeople);
%Generating the first set of assignments. The Gift Receivers
assignments = randperm(numPeople);


%This checks that no one is assigned to themselves and then regenerates
%this matrix if it
while any((numHolder(1:numPeople) == assignments(1:numPeople)) == 1 )
    assignments = randperm(numPeople);
end


%Sending email
for ii = 1:numPeople
    if ii ~= assignments(ii) && test == 0
        personName =  people{assignments(ii)}; %Grabs the gift receivers name.
        message = ['Ho ho ho, ' people{ii,1}, '  your Secret Santa is ', personName, '.']; %Creates the message.
        reminder = 'DO NOT LOSE THIS EMAIL - THIS IS THE ONLY COPY.';
        %Below sends the email.
        sendmail(people{ii,2}, 'Secret Santa Assignment', [message, ' ', reminder])
    elseif test == 1
       message = 'Tell Blank you got this';
       sendmail(people{ii,2}, 'Secret Santa Test', message)
    end
end
clear all %keeps it anonymous.    
