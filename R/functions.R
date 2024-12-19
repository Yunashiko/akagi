# データの読み込み ----

#' Load OTU Raw Data
#'
#' Reads the raw OTU data from a file and standardizes the column names.
#'
#' @param path The file path to the OTU data.
#' @return A tibble with cleaned OTU data.
#' @examples
#' load_otu_raw("data/otu_table.txt")
load_otu_raw <- function(path) {
  read_delim(
    path,
    delim = "\t", show_col_types = FALSE
  ) |>
    clean_names()
}

#' Load Plant Data
#'
#' Reads the plant data from a file and standardizes the column names.
#'
#' @param path The file path to the plant data.
#' @return A tibble with cleaned plant data.
#' @examples
#' load_plants("data/plants.csv")
load_plants <- function(path) {
  read_csv(path, show_col_types = FALSE) |>
    rename(
      sample = samplename,
      plant = Host.Plant,
      location = 採取地
    ) |>
    # サンプル名をOTUデータの標準化された名前に合わせる
    mutate(
      sample = make_clean_names(sample)
    )
}

# データクリーニング ----

#' Clean OTU Data
#'
#' Cleans the OTU data by removing low-frequency OTUs and standardizing the
#' data format.
#'
#' @param otu_raw A tibble with raw OTU data.
#' @param threshold_fraction The threshold fraction to filter low-frequency
#' OTUs.
#' @return A tibble with cleaned OTU data.
#' @examples
#' clean_otu_data(otu_raw, 0.001)
clean_otu_data <- function(otu_raw, threshold_fraction) {
  otu_raw |>
    # 以下の処理のため、データの形を縦長に変える
    pivot_longer(names_to = "host", values_to = "abun", cols = -otu_id) |>
    # 闘値をサンプルあたり総リードの0.1%と設定し、それより少ないリード数を0とみなす
    mutate(threshold = sum(abun) * threshold_fraction, .by = otu_id) |>
    mutate(abun = case_when(
      abun < threshold ~ 0,
      .default = abun
    )) |>
    select(-threshold) |>
    # 全ての植物において個体数の合計が0のOTUを除外する
    mutate(total_otu_abun = sum(abun), .by = otu_id) |>
    filter(total_otu_abun > 0) |>
    select(-total_otu_abun) |>
    # 縦長から幅広い形に戻す
    pivot_wider(names_from = host, values_from = abun)
}

#' Clean Taxonomy Data
#'
#' Cleans the taxonomy data by removing unnecessary characters and separating
#' the taxonomy levels into separate columns.
#'
#' @param otu_taxonomy A tibble with raw taxonomy data.
#' @return A tibble with cleaned taxonomy data.
#' @examples
#' clean_taxonomy_data(otu_taxonomy)
clean_taxonomy_data <- function(otu_taxonomy) {
  tax_levels <- c(
    "kingdom", "phylum", "class", "order", "family", "genus", "species"
  )
  otu_taxonomy |>
    mutate(taxonomy = str_remove_all(taxonomy, ".__")) |>
    separate_wider_delim(
      taxonomy,
      delim = "; ",
      names = tax_levels,
      too_few = "align_start"
    )
}

#' Count Taxonomy Ranks
#'
#' Counts the number of identified taxa at each taxonomic rank.
#'
#' @param taxonomy_clean A tibble with cleaned taxonomy data.
#' @return A tibble with the count of identified taxa at each rank.
#' @examples
#' count_tax_ranks(taxonomy_clean)
count_tax_ranks <- function(taxonomy_clean) {
  taxonomy_clean |>
    pivot_longer(
      names_to = "rank",
      values_to = "taxon",
      -otu_id
    ) |>
    filter(!is.na(taxon)) |>
    count(rank, sort = TRUE)
}

# 希釈 ----

#' Transpose OTU Data
#'
#' Transposes the OTU community matrix.
#'
#' @param otu_clean A tibble with cleaned OTU data.
#' @return A tibble with transposed OTU data.
#' @examples
#' transpose_otu(otu_clean)
transpose_otu <- function(otu_clean) {
  otu_clean |>
    pivot_longer(names_to = "plant", values_to = "abun", -otu_id) |>
    pivot_wider(names_from = "otu_id", values_from = "abun")
}

#' Untranspose OTU Data
#'
#' Untransposes the OTU community matrix.
#'
#' @param otu_clean A tibble with transposed OTU data.
#' @return A tibble with untransposed OTU data.
#' @examples
#' untranspose_otu(otu_clean)
untranspose_otu <- function(otu_clean) {
  otu_clean |>
    pivot_longer(names_to = "otu_id", values_to = "abun", -plant) |>
    pivot_wider(names_from = "plant", values_from = "abun")
}

#' Perform iNEXT Analysis
#'
#' Performs iNEXT analysis on the transposed OTU data.
#'
#' @param otu_clean_transposed A tibble with transposed OTU data.
#' @param nboot The number of bootstrap replicates.
#' @return The result of the iNEXT analysis.
#' @examples
#' i_next(otu_clean_transposed, nboot = 20)
i_next <- function(otu_clean_transposed, nboot = 20) {
  otu_clean_transposed |>
    untranspose_otu() |>
    column_to_rownames("otu_id") |>
    iNEXT(
      q = 0,
      nboot = nboot
    )
}

#' Extract Asymptotic Estimates
#'
#' Extracts asymptotic estimates from the iNEXT results.
#'
#' @param i_next_res A list of iNEXT results.
#' @return A tibble with cleaned asymptotic estimates.
#' @examples
#' extract_asy_est(i_next_res)
extract_asy_est <- function(i_next_res) {
  transpose(i_next_res) |>
    magrittr::extract2("AsyEst") |>
    bind_rows() |>
    as_tibble() |>
    janitor::clean_names()
}

#' Calculate Diversity Correlation
#'
#' Calculates the correlation between observed and estimated diversity.
#'
#' @param asymptotic_richness A tibble with asymptotic richness data.
#' @return A tibble with the correlation results.
#' @examples
#' calculate_diversity_corr(asymptotic_richness)
calculate_diversity_corr <- function(asymptotic_richness) {
  asymptotic_richness |>
    summarize(
      correlation = cor(observed, estimator),
      .by = diversity
    )
}

# ネットワーク解析 ----

#' Make Community Data for Network Analysis
#'
#' Prepares community data for network analysis by joining OTU, taxonomy, and
#' plant data.
#'
#' @param otu_clean A tibble with cleaned OTU data.
#' @param taxonomy_clean A tibble with cleaned taxonomy data.
#' @param plants A tibble with plant data.
#' @param tax_level The taxonomic level to use for the analysis.
#' @return A tibble with community data prepared for network analysis.
#' @examples
#' make_comm_for_network(otu_clean, taxonomy_clean, plants, "family")
make_comm_for_network <- function(
    otu_clean, taxonomy_clean, plants, tax_level) {
  otu_clean |>
    pivot_longer(names_to = "sample", values_to = "abun", -otu_id) |>
    left_join(taxonomy_clean, by = "otu_id", relationship = "many-to-one") |>
    left_join(plants, by = "sample", relationship = "many-to-one") |>
    filter(!is.na(.data[[tax_level]])) |>
    filter(abun > 0) |>
    mutate(abun = 1) |>
    select(all_of(tax_level), plant, abun) |>
    unique() |>
    pivot_wider(
      names_from = all_of(tax_level), values_from = abun, values_fill = 0
    )
}

#' Calculate Network Level Metrics
#'
#' Calculates network level metrics for the community data.
#'
#' @param network_comm A tibble with community data.
#' @param names_col The column name for the plant names.
#' @param ... Additional arguments passed to the networklevel function.
#' @return A list with network level metrics.
#' @examples
#' calc_network_level(network_comm)
calc_network_level <- function(network_comm, names_col = "plant", ...) {
  network_comm |>
    # 宿者が行、共生菌が列なので、このまま（転置行列化なしに）使う
    column_to_rownames(names_col) |>
    bipartite::networklevel(web = _, ...)
}

#' Randomize Single Community
#'
#' Randomizes a single community using the specified null model.
#'
#' @param network_comm A tibble with community data.
#' @param names_col The column name for the plant names.
#' @param null_model The null model to use for randomization.
#' @param n_iterations The number of iterations for randomization.
#' @return A randomized community.
#' @examples
#' randomize_single_comm(network_comm, null_model = "curveball", n_iterations = 1000)
randomize_single_comm <- function(
    network_comm,
    names_col = "plant", null_model = "curveball", n_iterations) {
  network_comm |>
    column_to_rownames(names_col) |>
    canaper::cpr_rand_comm(
      null_model = null_model, n_iterations = n_iterations
    )
}

#' Generate Random Communities
#'
#' Generates multiple randomized communities.
#'
#' @param network_comm A tibble with community data.
#' @param names_col The column name for the plant names.
#' @param null_model The null model to use for randomization.
#' @param n_iterations The number of iterations for randomization.
#' @param n_reps The number of replicates.
#' @return A list of randomized communities.
#' @examples
#' generate_rand_comms(network_comm, null_model = "curveball", n_iterations = 1000, n_reps = 10)
generate_rand_comms <- function(
    network_comm,
    names_col = "plant", null_model = "curveball", n_iterations, n_reps = 1) {
  replicate(
    n_reps,
    randomize_single_comm(
      network_comm,
      names_col = names_col,
      null_model = null_model,
      n_iterations = n_iterations
    ),
    simplify = FALSE
  )
}

#' Calculate P-Value
#'
#' Calculates the p-value for the observed value compared to random values.
#'
#' @param obs_val The observed value.
#' @param random_index_vals A vector of random index values.
#' @param index The index being tested.
#' @param species The species being tested (optional).
#' @return A tibble with the p-value and other statistics.
#' @examples
#' calculate_p_val(obs_val, random_index_vals, index)
calculate_p_val <- function(obs_val, random_index_vals, index, species = NULL) {
  tibble(
    obs = obs_val,
    # Calculate SES
    rand_mean = mean(random_index_vals, na.rm = TRUE),
    rand_sd = sd(random_index_vals, na.rm = TRUE),
    obs_z = (obs_val - rand_mean) / rand_sd,
    # Count number of times observed value is higher than random values
    obs_c_upper = canaper:::count_higher(obs_val, random_index_vals),
    # Count number of times observed value is lower than random values
    obs_c_lower = canaper:::count_lower(obs_val, random_index_vals),
    # Count the number of non-NA random values used for comparison
    obs_q = sum(!is.na(random_index_vals)),
    # Calculate p-value for upper tail
    obs_p_upper = obs_c_upper / obs_q,
    # Calculate p-value for lower tail
    obs_p_lower = obs_c_lower / obs_q
  ) |>
    mutate(species = species, .before = 0) |>
    mutate(index = index, .before = 0)
}

#' Run Randomization Test
#'
#' Runs a randomization test for the network community data.
#'
#' @param network_comm A tibble with community data.
#' @param random_comms A list of randomized communities.
#' @param names_col The column name for the plant names.
#' @param index The index being tested.
#' @return A tibble with the p-value and other statistics.
#' @examples
#' run_rand_test(network_comm, random_comms, index = "NODF")
run_rand_test <- function(
    network_comm,
    random_comms,
    names_col = "plant",
    index = "NODF") {
  obs_val <- calc_network_level(
    network_comm,
    names_col = names_col, index = index
  )

  random_index_vals <- purrr::map_dbl(
    random_comms, ~ networklevel(., index = index)
  )

  calculate_p_val(
    obs_val = obs_val, random_index_vals = random_index_vals,
    index = index
  )
}

#' Calculate D-Statistics
#'
#' Calculates D-statistics for the network community data.
#'
#' @param network_comm A tibble with community data.
#' @param names_col The column name for the plant names.
#' @param abuns Abundance data (optional).
#' @return A list with D-statistics.
#' @examples
#' calculate_d(network_comm)
calculate_d <- function(network_comm, names_col = "plant", abuns = NULL) {
  network_comm |>
    column_to_rownames(names_col) |>
    dfun(abuns = abuns)
}

#' Extract D-Statistics
#'
#' Extracts D-statistics from the random values.
#'
#' @param dfun_rand_vals A list of random D-statistics.
#' @param d_stat_select The D-statistic to select.
#' @return A list with the selected D-statistics.
#' @examples
#' extract_dstats(dfun_rand_vals, "dprime")
extract_dstats <- function(dfun_rand_vals, d_stat_select) {
  dfun_rand_vals |>
    magrittr::extract(
      str_detect(names(dfun_rand_vals), d_stat_select)) |>
    purrr::list_c()
}

#' Calculate P-Value for D-Statistics
#'
#' Calculates the p-value for the observed D-statistics compared to random
#' values.
#'
#' @param dfun_obs_res A list of observed D-statistics.
#' @param dfun_rand_vals A list of random D-statistics.
#' @param d_stat_select The D-statistic to select.
#' @param species_select The species to select.
#' @return A tibble with the p-value and other statistics.
#' @examples
#' calculate_p_val_dstats(dfun_obs_res, dfun_rand_vals, "dprime", "species1")
calculate_p_val_dstats <- function(
    dfun_obs_res, dfun_rand_vals,
    d_stat_select,
    species_select) {
  dfund_rand_vals_select <- extract_dstats(dfun_rand_vals, d_stat_select)

  calculate_p_val(
    obs_val = dfun_obs_res[[d_stat_select]][species_select],
    random_index_vals =
      dfund_rand_vals_select[names(dfund_rand_vals_select) == species_select],
    index = d_stat_select, species = species_select
  )
}

#' Check Binary Data
#'
#' Checks if some network indices always turn out 0 for binary data.
#'
#' @param index The network index to check.
#' @return The result of the network index calculation.
#' @examples
#' check_binary_data("weighted NODF")
check_binary_data <- function(index = "weighted NODF") {
  Safariland |>
    as.data.frame() |>
    tibble::rownames_to_column("plant") |>
    pivot_longer(names_to = "insect", values_to = "abun", -plant) |>
    as_tibble() |>
    mutate(abun = if_else(abun > 0, 1, 0)) |>
    pivot_wider(names_from = insect, values_from = abun) |>
    column_to_rownames("plant") |>
    bipartite::networklevel(index = index)
}

#' Make Network Graph
#'
#' Creates a network graph from the community data.
#'
#' @param network_comm A tibble with community data.
#' @return A tbl_graph object representing the network graph.
#' @examples
#' make_network_graph(network_comm)
make_network_graph <- function(network_comm) {
  network_comm |>
    column_to_rownames("plant") |>
    as.matrix() |>
    graph_from_biadjacency_matrix(multiple = TRUE) |>
    as_tbl_graph() |>
    activate(nodes) |>
    # 図で示すデータを用意する
    # - 植物かどうか
    rename(plant = type) |>
    mutate(plant = !plant) |>
    # - 媒介中心性
    mutate(importance = centrality_betweenness()) |>
    # - 植物の名前
    mutate(plant_name = if_else(plant, name, NA_character_))
}