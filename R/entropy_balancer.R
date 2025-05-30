#' @title Run Entropy Balancer for Multiple Zones
#' @description This function runs two-stage entropy balancing on survey seed data. It performs:
#'   (1) initial zone-level balancing,
#'   (2) meta control factoring from the results,
#'   (3) final balancing with the meta targets.
#' @param seed_df A data.frame containing the seed incidence and weight data.
#' @param controls_zones_df Data.frame of per-zone control totals (zone_id, control columns).
#' @param controls_region_df Optional. Data.frame of per-zone meta control totals.
#' @param controls_importance Optional. Vector of control importances.
#' @param sample_id_col The name of the sample ID column.
#' @param zone_id_col The name of the zone ID column.
#' @param region_id_col The name of the region ID column.
#' @param weight_col The name of the initial weight column.
#' @param lb Lower bound vector or scalar factor for weights (default = 0.1).
#' @param ub Upper bound vector or scalar factor for weights (default = 4).
#' @param lb_abs Absolute lower bound for weights (default = 0).
#' @param ub_abs Absolute upper bound for weights (default = Inf).
#' @param max_iterations Maximum number of iterations (default = 10000).
#' @param max_delta Convergence threshold (default = 1e-8).
#' @param master_control Index of master control, or -1 for none.
#' @param print_every_n Iteration progress frequency.
#' @return A list with round1, meta_controls, and round2 results.
#' @export
#' @import data.table
#' @import Rcpp
#' @useDynLib rentropy, .registration = TRUE
multi_entropy_balancer = function(
    seed_df,
    controls_zones_df,
    controls_region_df = NULL,
    controls_importance = NULL,    
    sample_id_col = "sample_id",
    zone_id_col = "zone_id",
    region_id_col = "region_id",
    weight_col = "initial_weight",
    lb = 0.1,
    ub = 4,
    lb_abs = 0,
    ub_abs = Inf,
    max_iterations = 10000,
    max_delta = NULL,
    master_control = NULL,
    print_every_n = 100
) {

    # Convert to data.table
    seed_df = as.data.table(seed_df)
    controls_zones_df = as.data.table(controls_zones_df)
    zone_ids = unique(seed_df[[zone_id_col]])
    control_names = setdiff(colnames(seed_df), c(sample_id_col, zone_id_col, region_id_col, weight_col))
    master_control_index = if (is.null(master_control)) -1 else match(master_control, control_names)


    # If Null, aggregate controls_zones_df to per-zone totals
    if (!is.null(controls_region_df)) {
        controls_region = as.matrix(controls_region_df)[1, control_names, drop = FALSE]
    } else {
        controls_region = colSums(controls_zones_df[, control_names, with = FALSE], na.rm = TRUE)
    }

    # Build per-zone data
    id_list = list()
    incidence_list = list()
    initial_weights_list = list()
    control_zones_list = list()
    lb_list = list()
    ub_list = list()

    # Prepare list data
    for (zone in zone_ids) {
        sub_seed = seed_df[get(zone_id_col) == zone]
        sub_control = controls_zones_df[get(zone_id_col) == zone]

        incidence_list[[zone]] = as.matrix(sub_seed[, control_names, with = FALSE])
        initial_weights_list[[zone]] = sub_seed[[weight_col]]
        control_zones_list[[zone]] = as.numeric(sub_control[, control_names, with = FALSE])
        id_list[[zone]] = sub_seed[, c(sample_id_col, zone_id_col, region_id_col), with = FALSE]

        # Prepare weight bounds
        if (length(lb) == 1) {
            lb_ = initial_weights_list[[zone]] * lb
        } else if (length(lb) != length(initial_weights_list[[zone]])) {
            stop("Length of 'lb' must match number of controls.")
        }


        if (length(ub) == 1) {
            ub_ = initial_weights_list[[zone]] * ub
        } else if (length(ub) != length(initial_weights_list[[zone]])) {
            stop("Length of 'ub' must match number of controls.")
        }

        # Cap the weight bounds. Default is 0 and Inf, use whichever is narrower.
        lb_list[[zone]] = pmax(lb_, lb_abs)
        ub_list[[zone]] = pmin(ub_, ub_abs)
    }

    # Apply importance or default
    if (is.null(controls_importance)) {
        importance = rep(1.0, length(control_names))

    } else if (length(controls_importance) == 1 && is.numeric(controls_importance)) {
        importance = rep(controls_importance, length(control_names))
    
    } else if (length(controls_importance) == length(control_names)) {
        if (is.null(names(controls_importance)) || any(is.na(names(controls_importance)))) {
            stop("If 'controls_importance' has length > 1, it must be a named vector.")
        }
        importance = as.numeric(controls_importance[control_names])

    } else if (length(controls_importance) != length(control_names)) {
        stop("Length of 'importance' must match number of controls.")

    }


    # --- ROUND 1 ---
    results = lapply(seq_along(zone_ids), function(i) {
        print(paste0("Processing zone: ", i))
        res = entropy_balancer(
            incidence            = incidence_list[[i]],
            initial_weights      = initial_weights_list[[i]],
            control_totals       = control_zones_list[[i]],
            controls_importance  = importance,
            weights_lb           = lb_list[[i]],
            weights_ub           = ub_list[[i]],
            master_control_index = master_control_index,
            max_iterations       = max_iterations,
            max_delta            = max_delta,
            print_every_n        = print_every_n
        )
    })

    # Combine the results with ID list
    final_weights = rbindlist(
        lapply(seq_along(results), function(i) {
            cbind(id_list[[i]], weight = results[[i]]$weights)
        })
    )

    return(final_weights)
}


# NOT USED?
#' @title Meta Control Factoring (wide control tables)
#' @description Computes zone-level meta controls by scaling weighted incidence to match region-level targets.
#' @param seed_weights_dt Data.table with seed incidence and weights.
#' @param meta_controls_df Data.frame of per-region meta control totals (region_id, control columns).
#' @param zone_id_col Name of the zone ID column in seed_weights_dt.
#' @param region_id_col Name of the region ID column in seed_weights_dt.
#' @param weight_col Name of the weight column in seed_weights_dt.
#' @return A data.table with zone_id and scaled meta controls.
meta_control_factoring = function(
    seed_weights_dt,
    meta_controls_df,
    zone_id_col,
    region_id_col,
    weight_col
) {
    meta_controls = names(meta_controls_df)

    # Apply weight to seed
    seed_weights_dt_ = copy(seed_weights_dt)
    for (target in meta_controls) {
        seed_weights_dt_[, (target) := get(target) * get(weight_col)]
    }

    # Sum by zone
    zone_totals = seed_weights_dt_[
        ,
        lapply(.SD, sum, na.rm = TRUE),
        by = zone_id_col,
        .SDcols = meta_controls
    ]

    # Compute scaling factors
    total_controls = seed_weights_dt_[, lapply(.SD, sum, na.rm = TRUE), .SDcols = meta_controls]

    # Apply scaling and round
    scale_factors = meta_controls_df / total_controls

    # Replace NA with 1
    scale_factors = scale_factors[, lapply(.SD, function(x) ifelse(!is.finite(x), 1, x))]

    # Apply scaling factors to zone totals
    for (target in meta_controls) {
        zone_totals[, (target) := round(get(target) * scale_factors[[target]])]
    }

    # Coerce to integer
    for (target in meta_controls) {
        set(zone_totals, j = target, value = as.integer(zone_totals[[target]]))
    }

    # Split data frame into a list of row vectors    
    zone_totals_list = split(zone_totals, zone_totals[[zone_id_col]])
    for (i in seq_along(zone_totals_list)) {
        # Drop zone_id_col and convert to matrix
        dt = zone_totals_list[[i]][, !zone_id_col, with = FALSE]
        zone_totals_list[[i]] = as.matrix(dt)[1,]
    }

    return(zone_totals_list)
}
