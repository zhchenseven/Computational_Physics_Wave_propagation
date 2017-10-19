from Make_movie import Gen_movie
import matplotlib.pyplot as plt
import numpy as np

Plot=Gen_movie()
M=30
fd_name='shock_DN_F1_UW1'
index=0
data_type='dat'
size_sufx='size'
main_fd_name='shock_DN_F1_scan_T_0p5'

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
title_str='shock Dirichlet and Neumann F1 flux when aT=0.5'
err_arr = np.array([[0.1756998412,0.1159573302,0.07437298647,0.02941176471],
                    [0.07926643738,0.062147274,0.04552373933,0.02941176471], \
                    [0.1065563774,0.0679584965,0.04310658066,0.02941176471],
                    [0.08325454092,0.0679584965,0.04783390642,0.02941176471]])
Plot.plot_err_norm(nu,err_arr,title_str,scheme_str)


plt.show()