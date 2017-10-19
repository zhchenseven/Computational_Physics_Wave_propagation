from Make_movie import Gen_movie


Movie=Gen_movie()
M=17
fd_name='shock_DN_F1_scan_T_0p5_0/shock_DN_F1_UW1_nu_0p5'
index=0
data_type='dat'
size_sufx='size'

Movie.set_basic(M,fd_name,data_type,size_sufx,index)

x_name='R1_x'
t_out_name='R1_t_out'
u_exa_name='R2_u_exact'
u_num_name='R2_u_track'

mov_name='demo'
fps=5
fine_fac=4
Movie.set_movie_para(mov_name,fps,fine_fac)

Movie.set_var(x_name,t_out_name,u_num_name,u_exa_name)
Movie.read_raw_data()
Movie.gen_movie()