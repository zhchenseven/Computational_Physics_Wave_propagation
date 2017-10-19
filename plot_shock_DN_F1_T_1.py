from Make_movie import Gen_movie
import matplotlib.pyplot as plt
import numpy as np

Plot=Gen_movie()
M=30
fd_name='shock_DN_F1_UW1'
index=0
data_type='dat'
size_sufx='size'
main_fd_name='shock_DN_F1_scan_T_1'

Plot.set_basic(M,fd_name,data_type,size_sufx,index)
main_fd_ctrl=1
Plot.set_main_fd_ctrl(main_fd_ctrl,main_fd_name)

x_name='R1_x'
t_out_name='R1_t_out'
u_exa_name='R2_u_exact'
u_num_name='R1_uf'

mov_name='demo'
fps=10
fine_fac=4
Plot.set_movie_para(mov_name,fps,fine_fac)

Plot.set_var(x_name,t_out_name,u_num_name,u_exa_name)
# Plot.read_raw_data()
# Plot.gen_movie()
scheme_str=['LF','LW','MD','UW1']
fd_name_base_o='shock_DN_F1'
var_name='nu'
var_str=['1','0p75','0p5','0p25']
var_val=[1,0.75,0.5,0.25]
y_min=-0.1
y_max=1.3
for i in range(len(scheme_str)):
    fd_name_base=fd_name_base_o+'_'+scheme_str[i]
    Plot.plot_group_data(fd_name_base,var_name,var_str,var_val,y_min,y_max)

nu=[0.25,0.5,0.75,1]
title_str='shock Dirichlet and Neumann F1 flux when aT=1'
err_arr = np.array([[0.04641024655,0.01648920167,0.003602943027,1e-10],
                    [0.02066685788,0.008694746278,0.003020102171,1e-10], \
                    [0.01169086506,0.001564054043,5.464146928e-06,1e-10],
                    [0.004017098275,0.001564054043,0.0001886489281,1e-10]])
log_scale_ctrl=1
Plot.plot_err_norm(nu,err_arr,title_str,scheme_str,1e-11,1,1)


plt.show()