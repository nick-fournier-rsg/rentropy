# testthat::test_dir("tests/testthat", filter = "balancer_rcpp")

# To build the Rcpp code, run:
# Rcpp::compileAttributes(); devtools::document()

testthat::test_that("entropy_balancer handles large arbitrary-sized cases efficiently", {
    # A typical survey size is around 5k. A large one maybe 15k.
    set.seed(123)
    n_rows = 2500
    n_cols = 50

    # Generate a random incidence table (0-2 values)
    incidence_table = matrix(
        sample(0:2, n_rows * n_cols, replace = TRUE),
        nrow = n_rows, ncol = n_cols
    )
    initial_weights = rep(1.0, n_rows)
    controls = colSums(incidence_table) * runif(n_cols, 0.8, 1.2)
    importance = rep(1.0, n_cols)
    weights_lb = rep(0.1, n_rows)
    weights_ub = rep(30.0, n_rows)

    # Time the balancer
    result = entropy_balancer(
        incidence = incidence_table,
        initial_weights = initial_weights,
        control_totals = controls,
        controls_importance = importance,
        weights_lb = weights_lb,
        weights_ub = weights_ub,
        master_control_index = -1,
        max_iterations = 10000,
        max_delta = 1e-8,
        print_every_n = 50
    )

    # Check that result is a list and contains weights
    testthat::expect_type(result, "list")
    testthat::expect_true("weights" %in% names(result))

    # Check that weights are numeric and of correct length
    testthat::expect_type(result$weights, "double")
    testthat::expect_length(result$weights, n_rows)

    # Check that weighted sums are close to controls
    weighted_sum = unname(round(colSums(incidence_table * result$weights), 2))
    testthat::expect_equal(weighted_sum, round(controls, 2), tolerance = 0.1)

    # Check that runtime is reasonable (e.g., < 10 seconds)
    testthat::expect_lt(result$elapsed_time, 10)
})

testthat::test_that("single zone balancer for simple case", {

    incidence_table = as.matrix(data.frame(
        hh_1 = c(1, 1, 1, 0, 0, 0, 0, 0),
        hh_2 = c(0, 0, 0, 1, 1, 1, 1, 1),
        p1   = c(1, 1, 2, 1, 0, 1, 2, 1),
        p2   = c(1, 0, 1, 0, 2, 1, 1, 1),
        p3   = c(1, 1, 0, 2, 1, 0, 2, 0)
    ))

    # one weight per row in incidence table
    initial_weights = rep(1.0, nrow(incidence_table))

    # column totals which the final weighted incidence table sums must satisfy
    controls = c(35, 65, 91, 65, 104)

    # one for every column in incidence_table
    importance = rep(100000.0, ncol(incidence_table))
    
    result = entropy_balancer(
        incidence = incidence_table,
        initial_weights = initial_weights,
        control_totals = controls,
        controls_importance = importance,
        weights_lb = rep(0.0, length(initial_weights)),
        weights_ub = rep(30.0, length(initial_weights)),
        master_control_index = -1,
        max_iterations = 10000,
        max_delta = 1e-8,
        print_every_n = 10
    )

    expected_weights = c(8.939684, 23.443184, 2.615558, 25.898310, 14.346718, 11.007512, 2.736926, 11.007512)

    published_final_weights = c(1.36, 25.66, 7.98, 27.79, 18.45, 8.64, 1.47, 8.64)
    published_weighted_sum = unname(round(colSums(incidence_table * published_final_weights), 2))
    weighted_sum = unname(round(colSums(incidence_table * result$weights), 2))

    testthat::expect_equal(result$weights, expected_weights, tolerance = 1e-4)
    testthat::expect_equal(weighted_sum, published_weighted_sum, tolerance = 0.1)
    testthat::expect_equal(weighted_sum, controls, tolerance = 0.1)
})

# Rcpp::compileAttributes(); devtools::document()
testthat::test_that("Check results against populationsim", {

    # Read in tables
    # Get fixture from inst/fixtures
    fixture_path = system.file("test_fixtures", package = "rentropy")

    seed_df = fread(file.path(fixture_path, "incidence_table.csv"), header = TRUE)
    controls_zones_df = fread(file.path(fixture_path, "controls.csv"), header = TRUE)
    importance = fread(file.path(fixture_path, "importance.csv"), header = TRUE)
    controls_importance = setNames(importance$importance, importance$target)
    expected_weights = fread(file.path(fixture_path, "weights.csv"), header = TRUE)

    result = multi_entropy_balancer(
        seed_df = seed_df,
        controls_zones_df = controls_zones_df,
        controls_importance = controls_importance,
        sample_id_col = "hh_id",
        zone_id_col = "SUBREGCluster",
        region_id_col = "Region",
        weight_col = "sample_weight",
        lb = 0.5,
        ub = 4.0,
        lb_abs = 1,
        ub_abs = 20000,
        master_control = "HH_Total",
        max_iterations = 10000,
        print_every_n = 500,
        max_delta = 1e-8
    )

    # Check the results
    testthat::expect_equal(
        expected_weights[order(hh_id), weight],
        result[order(hh_id), weight],
        tolerance = 0.1
    )

    # Plot check
    # plot(
    #     x = expected_weights[order(hh_id), weight],
    #     y = result[order(hh_id), weight],
    # )

    # Check that the result is close to the target totals
    testthat::expect_equal(sum(result$weight), sum(controls_zones_df$HH_Total))

    # Check that result is a list and contains weights
    testthat::expect_type(result, "list")
    testthat::expect_true("weight" %in% names(result))

    # Check that weights are numeric and of correct length
    testthat::expect_type(result$weight, "double")
    testthat::expect_length(result$weight, nrow(seed_df))

    # Get sum by zone
    weight_zones_sums = cbind(
        seed_df[, .(hh_id, SUBREGCluster)],
        seed_df[, names(controls_importance), with = FALSE] * result$weight
    )
    weight_zones_sums = weighted_zones_sums[,
        lapply(.SD, sum),
        .SDcols = names(controls_importance),
        by = SUBREGCluster
    ]

    # Check that weighted sums are close to controls
    y_ = weight_zones_sums[, names(controls_importance), with = FALSE]
    y = controls_zones_df[, names(controls_importance), with = FALSE]
    n = ncol(y_) * nrow(y_)
    avg = mean(as.matrix(y))

    # Pretty bad, but so was the original populationsim result
    sqrt(sum((y_ - y)^2) / n) / avg
})
