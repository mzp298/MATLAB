function mail2me(mail,subject,content,DataPath)
try  %skip when no internet
MailAddress = 'mzp298@gmail.com';
password = 'Mzp20177';
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.gmail.com'); 
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail(mail,subject,content,DataPath);
catch
end
end
% mail2me('mzp298@gmail.com','Job done','b=1.02','F:\Git\MATLAB\anew\alpha_varying.xlsx');