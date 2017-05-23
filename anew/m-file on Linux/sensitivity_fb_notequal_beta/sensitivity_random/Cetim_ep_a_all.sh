#!/bin/bash
screen -dmS a_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_a_max.m" 
#!/bin/bash
screen -dmS a_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_a_min.m"
#!/bin/bash
screen -dmS beta_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_beta_max.m"
#!/bin/bash
screen -dmS beta_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_beta_min.m"
#!/bin/bash
screen -dmS lam_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_lam_max.m"
#!/bin/bash
screen -dmS lam_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_lam_min.m"
#!/bin/bash
screen -dmS W0_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_W0_max.m"
#!/bin/bash
screen -dmS W0_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_W0_min.m"
#!/bin/bash
screen -dmS pb_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_pb_max.m"
#!/bin/bash
screen -dmS pb_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_notequal_beta/sensitivity_random/Cetim_ep_a_05_pb_min.m"



