source('0_helpers.R')
source('1_simulations.R')
source('2_model_comparison.R')
source('3_hat_matrix.R')
source('4_principaled_mse.R')
source('5_plots.R')

#==================================================#
#   GENERATING AND SAVING RESULTS                  #
#==================================================#
set.seed(64)

# generate and save the posterior draws to /posterior
fn_save_post()

# generate and save the hat matrix weights results to /posterior
fn_save_hat_matrix()

#==================================================#
#   TABLES                                         #
#==================================================#
set.seed(867)

# results for model comparison
model_comp<-fn_model_comparison()

# output for Table 3
fn_table3(model_comp)

# output for Table 4
fn_table4(model_comp)

# output for Table 5
fn_table5()

#==================================================#
#   FIGURES                                        #
#==================================================#
set.seed(5309)

measures<-fn_data_prep()
post<-fn_load_post(measures)

# output for Figure 1
fn_figure1()

# output for Figure 2
fn_figure2(measures,post)

# output for Figure 3
fn_figure3(measures,post)

# output for Figure 4
fn_figure4(measures,post)
