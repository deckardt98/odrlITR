test_that("tree depth rejects invalid constructor values", {
  expect_error(odrl_control(tree_depth = -1L), "`tree_depth`")
  expect_error(odrl_control(tree_depth = 1.5), "`tree_depth`")
  expect_error(odrl_control(tree_depth = c(1L, 2L)), "`tree_depth`")
  expect_error(odrl_control(tree_depth = NA_integer_), "`tree_depth`")

  expect_identical(odrl_control(tree_depth = 0L)$tree_depth, 0L)
  expect_identical(odrl_control(tree_depth = 3L)$tree_depth, 3L)
})

test_that("policytree options cannot override dedicated tree controls", {
  skip_if_not_installed("policytree")

  x <- matrix(c(-2, -1, 1, 2), ncol = 1)
  score <- c(-2, -1, 1, 2)
  reserved <- c("X", "Gamma", "depth", "split.step", "min.node.size")

  for (option in reserved) {
    control <- odrl_control(
      tree_depth = 1L,
      tree_min_node_size = 1L,
      tree_options = stats::setNames(list(1), option)
    )
    expect_error(
      odrlITR:::.odrl_fit_tree(x, score, control),
      "cannot override dedicated control field",
      fixed = TRUE,
      info = option
    )
  }
})

test_that("custom tree diagnostics persist resolved search controls", {
  backend <- list(
    name = "diagnostic_tree",
    fit = function(x, score, rewards, depth, min_node_size, split_step,
                   options) {
      list(
        depth = depth,
        min_node_size = min_node_size,
        split_step = split_step,
        options = options
      )
    },
    predict = function(model, newx) ifelse(newx[, 1L] >= 0, 1, -1)
  )
  control <- odrl_control(
    tree_depth = 3L,
    tree_min_node_size = 7L,
    tree_split_step = 4L,
    tree_backend = backend,
    tree_options = list(label = "audit")
  )
  x <- matrix(c(-2, -1, 1, 2), ncol = 1)
  fit <- odrlITR:::.odrl_fit_tree(x, c(-2, -1, 1, 2), control)

  expect_identical(fit$engine, "diagnostic_tree")
  expect_identical(fit$fit$depth, 3L)
  expect_identical(fit$fit$min_node_size, 7L)
  expect_identical(fit$fit$split_step, 4L)
  expect_identical(fit$fit$options, list(label = "audit"))
  expect_identical(fit$diagnostics$depth, 3L)
  expect_identical(fit$diagnostics$min_node_size, 7L)
  expect_identical(fit$diagnostics$split_step, 4L)
  expect_identical(fit$diagnostics$options, list(label = "audit"))
  expect_identical(fit$diagnostics$backend, "diagnostic_tree")
  expect_false(fit$diagnostics$global_candidate_search)
  expect_match(fit$diagnostics$search_scope, "depth <= 3", fixed = TRUE)
  expect_match(fit$diagnostics$search_scope, "minimum node size 7", fixed = TRUE)
  expect_match(fit$diagnostics$search_scope, "split.step 4", fixed = TRUE)
})
