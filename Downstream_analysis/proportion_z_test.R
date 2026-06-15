obs_AA = 88
obs_BB = 58
obs_AB = 37
obs_BA = 37
total_obs = obs_AA + obs_BB + obs_AB + obs_BA

exp_AA = 36408
exp_BB = 151548
exp_AB = 74402
exp_BA = 74402
total_exp = exp_AA + exp_BB + exp_AB + exp_BA

obs_homo_prop  <- (obs_AA + obs_BB) / total_obs
obs_heter_prop <- (obs_AB + obs_BA) / total_obs

exp_homo_prop  <- (exp_AA + exp_BB) / total_exp
exp_heter_prop <- (exp_AB + exp_BA) / total_exp

prop.test(
  x = c(obs_AA + obs_BB, exp_AA + exp_BB),
  n = c(total_obs, total_exp)
)


obs_AA = 28
obs_BB = 18
obs_AB = 8
obs_BA = 8
total_obs = obs_AA + obs_BB + obs_AB + obs_BA

exp_AA = 42182
exp_BB = 140476
exp_AB = 77051
exp_BA = 77051
total_exp = exp_AA + exp_BB + exp_AB + exp_BA

obs_homo_prop  <- (obs_AA + obs_BB) / total_obs
obs_heter_prop <- (obs_AB + obs_BA) / total_obs

exp_homo_prop  <- (exp_AA + exp_BB) / total_exp
exp_heter_prop <- (exp_AB + exp_BA) / total_exp

prop.test(
  x = c(obs_AA + obs_BB, exp_AA + exp_BB),
  n = c(total_obs, total_exp)
)

