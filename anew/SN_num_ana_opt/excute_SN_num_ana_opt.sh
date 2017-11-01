#!/bin/bash
screen -dmS SN_num_ana_30steps_J2 matlab7 -nosplash -nodesktop -r "run /home/ma/SN_num_ana_opt/SN_num_ana_30steps_J2.m" 
#!/bin/bash
screen -dmS SN_num_ana_100steps_J2 matlab7 -nosplash -nodesktop -r "run /home/ma/SN_num_ana_opt/SN_num_ana_100steps_J2.m" 
#!/bin/bash
screen -dmS SN_opt_ana_dalp000001_200steps_J2 matlab7 -nosplash -nodesktop -r "run /home/ma/SN_num_ana_opt/SN_opt_ana_dalp000001_200steps_J2.m" 
#!/bin/bash
screen -dmS SN_opt_ana_dalp000002_200steps_J2 matlab7 -nosplash -nodesktop -r "run /home/ma/SN_num_ana_opt/SN_opt_ana_dalp000002_200steps_J2.m" 
