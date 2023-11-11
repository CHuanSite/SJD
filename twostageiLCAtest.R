dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
comp_num = c(2,2,2,2,2,2,2,2,2)
proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
res_twoStageiLCA = twoStageiLCA(
dataset,
group,
comp_num,
proj_dataset = proj_dataset,
proj_group = proj_group)



weighting = NULL
backup = 0
screen_prob = NULL
enable_normalization = TRUE
column_sum_normalization = FALSE