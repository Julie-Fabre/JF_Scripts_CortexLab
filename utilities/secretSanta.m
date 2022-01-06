
%For this to work in Gmail you must allow less secure connections in gmail.
setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
setpref('Internet', 'E_mail', 'username@gmail.com');
setpref('Internet', 'SMTP_Username', 'username');
setpref('Internet', 'SMTP_Password', 'password');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth', 'true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port', '465');

people = {};
people{1, 1} = 'Julie';
people{1, 2} = 'fakeemail4@gmail.com';

numPeople = size(people, 1);
numHolder = linspace(1, numPeople, numPeople);
assignments = randperm(numPeople);

%This checks that no one is assigned to themselves and then regenerates
%this matrix if it
while any((numHolder(1:numPeople) == assignments(1:numPeople)) == 1)
    assignments = randperm(numPeople);
end

%Sending email
for iPerson = 1:numPeople

    personName = people{assignments(iPerson)};
    message = ['Ho ho ho, ', people{iPerson, 1}, '  your Secret Santa is ', personName, '.'];
    reminder = 'DO NOT LOSE THIS EMAIL - THIS IS THE ONLY COPY.';
    sendmail(people{iPerson, 2}, 'Secret Santa Assignment', [message, ' ', reminder])

end
clearvars assignments; %keeps it anonymous.
