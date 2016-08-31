function mail2me(subject,content) 
MailAddress ='mzp298@gmail.com';
password = 'Mzp20099';
setpref('Internet','E_mail',MailAddress); 
setpref('Internet','SMTP_Server','smtp.gmail.com'); 
setpref('Internet','SMTP_Username',MailAddress); 
setpref('Internet','SMTP_Password',password); 
props = java.lang.System.getProperties; 
props.setProperty('mail.smtp.auth','true'); 
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail('mzp298@gmail.com',subject,content);