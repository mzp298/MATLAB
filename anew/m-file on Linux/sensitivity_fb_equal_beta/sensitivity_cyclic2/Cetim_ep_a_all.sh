#!/bin/bash
screen -dmS a_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_a_max.m" 
#!/bin/bash
screen -dmS a_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_a_min.m"
#!/bin/bash
screen -dmS beta_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_beta_max.m"
#!/bin/bash
screen -dmS beta_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_beta_min.m"
#!/bin/bash
screen -dmS lam_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_lam_max.m"
#!/bin/bash
screen -dmS lam_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_lam_min.m"
#!/bin/bash
screen -dmS W0_max matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_W0_max.m"
#!/bin/bash
screen -dmS W0_min matlab -nosplash -nodesktop -r "run /home/ma/sensitivity_fb_equal_beta/sensitivity_cyclic2/Cetim_ep_a_01_W0_min.m"


