from Make_movie import Gen_movie
import matplotlib.pyplot as plt
import numpy as np

Plot=Gen_movie()
M=30
fd_name='sin_PD_F1_UW1'
index=0
data_type='dat'
size_sufx='size'
main_fd_name='cos_PD_F2_scan'

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
fd_name_base_o='cos_PD_F2'
var_name='nu'
var_str=['1','0p75','0p5','0p25']
var_val=[1,0.75,0.5,0.25]
y_min=-0.2
y_max=1.2
for i in range(len(scheme_str)):
    fd_name_base=fd_name_base_o+'_'+scheme_str[i]
    Plot.plot_group_data(fd_name_base,var_name,var_str,var_val,y_min,y_max)

nu = [0.25, 0.5, 0.75, 1]
title_str = 'cosine periodic F2 flux'
err_arr = np.array([[0.05815778044,0.05366503119,0.04596604384,0.03252356717]])
scheme_str=['Flux 2']
Plot.plot_err_norm(nu, err_arr, title_str, scheme_str)

ampd_err = np.array([[0.724997835,0.7664966715,0.8384281541,1.007122486]])
phase_err = np.array([[0.0, 0.030303030303030276, 0.030303030303030276, 0.030303030303030276]])

phase_err*=2*np.pi
err_type='ampd'
dx=0.0303030303
km=np.pi/(2*dx)
theta=km*dx
plot_ana_ctrl=0
Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,ampd_err,err_type,plot_ana_ctrl)
err_type='phase'
Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,phase_err,err_type,plot_ana_ctrl,-0.015,0.07)


plt.show()