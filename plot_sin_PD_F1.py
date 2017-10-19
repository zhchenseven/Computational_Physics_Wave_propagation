from Make_movie import Gen_movie
import matplotlib.pyplot as plt
import numpy as np

Plot=Gen_movie()
M=30
fd_name='sin_PD_F1_UW1'
index=0
data_type='dat'
size_sufx='size'
main_fd_name='sin_PD_F1_scan'

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
fd_name_base_o='sin_PD_F1'
var_name='nu'
var_str=['1','0p75','0p5','0p25']
var_val=[1,0.75,0.5,0.25]
y_min=-1.1
y_max=1.1
for i in range(len(scheme_str)):
    fd_name_base=fd_name_base_o+'_'+scheme_str[i]
    Plot.plot_group_data(fd_name_base,var_name,var_str,var_val,y_min,y_max)

# now plot the L1 norm error
nu=[0.25,0.5,0.75,1]
title_str='sinusoidal periodic F1 flux'
err_arr=np.array([[0.5967723905,0.4589340128,0.2529912857,2.016360934e-16],[0.03444765807,0.02755175037,0.01606591351,2.016360934e-16],\
                  [0.4163716944,0.2234622603,0.09915404499,2.016360934e-16],[0.3032028344,0.2234622603,0.1244506931,2.016360934e-16]])
ampd_err=np.array([[0.03398046218,0.2598977031,0.5944863376,1],[0.9952356331,0.9940539705,0.9958724535,1]\
                      ,[0.3256369889,0.6380761394,0.8394126754,1],[0.5094467491,0.6380761394,0.7992518172,1]])
# phase_err=np.array([[0.7575757576,0.7575757576,0.7575757576,0.7272727273],[0.7272727273,0.7575757576,0.7575757576,0.7272727273]\
#                        ,[0.7575757576,0.7575757576,0.7575757576,0.7272727273],[0.7575757576,0.7575757576,0.7575757576,0.7272727273]])
phase_err=np.array([[0.517131,0.513711,0.50793,0.484848],[0.491542,0.493249,0.496076,0.484848]\
                       ,[0.499999,0.5,0.500006,0.484848],[0.496603,0.5,0.511106,0.484848]])
N_itr=np.array([198,99,66,49])
dx=0.0303030303
phase_err*=2*np.pi

for i in range(4):
    for j in range(4):
        phase_err[i,j]-=phase_err[i,3]
# phase_err-=dx*np.pi
Plot.plot_err_norm(nu,err_arr,title_str,scheme_str,-0.03,0.65)
err_type='ampd'

km=2*np.pi
theta=km*dx
plot_ana_ctrl=1
# Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,ampd_err,err_type,plot_ana_ctrl,0,1.05)
Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,ampd_err,err_type,plot_ana_ctrl,0,1.1,N_itr)
err_type='phase'
scheme_str=['LF','LW','MD','UW1']
# Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,phase_err,err_type,plot_ana_ctrl,-0.002,0.012)
Plot.plot_ampd_phase_err(nu,theta,scheme_str,title_str,phase_err,err_type,plot_ana_ctrl,-0.05,0.1,N_itr)
plt.show()