# データの読み込み ----

# 生データを読み込んで、列名（サンプル名）を標準化する
load_otu_raw <- function(path) {
  read_delim(
    path,
    delim = "\t", show_col_types = FALSE
  ) |>
    clean_names()
}

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

# Transpose the OTU community matrix
# 群集マトリックスを転置行列化する
#
# Normally columns are OTU and rows are plant.
# This reverses that.
transpose_otu <- function(otu_clean) {
  otu_clean |>
    pivot_longer(names_to = "plant", values_to = "abun", -otu_id) |>
    pivot_wider(names_from = "otu_id", values_from = "abun")
}

untranspose_otu <- function(otu_clean) {
  otu_clean |>
    pivot_longer(names_to = "otu_id", values_to = "abun", -plant) |>
    pivot_wider(names_from = "plant", values_from = "abun")
}

i_next <- function(otu_clean_transposed, nboot = 20) {
  otu_clean_transposed |>
    untranspose_otu() |>
    column_to_rownames("otu_id") |>
    iNEXT(
      q = 0,
      nboot = nboot
    )
}

extract_asy_est <- function(i_next_res) {
  transpose(i_next_res) |>
    magrittr::extract2("AsyEst") |>
    bind_rows() |>
    as_tibble() |>
    janitor::clean_names()
}

calculate_diversity_corr <- function(asymptotic_richness) {
  asymptotic_richness |>
    summarize(
      correlation = cor(observed, estimator),
      .by = diversity
    )
}

# ネットワーク解析 ----

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

calc_network_level <- function(network_comm, names_col = "plant", ...) {
  network_comm |>
    # 宿者が行、共生菌が列なので、このまま（転置行列化なしに）使う
    column_to_rownames(names_col) |>
    bipartite::networklevel(web = _, ...)
}

randomize_single_comm <- function(
    network_comm,
    names_col = "plant", null_model = "curveball", n_iterations) {
  network_comm |>
    column_to_rownames(names_col) |>
    canaper::cpr_rand_comm(
      null_model = null_model, n_iterations = n_iterations
    )
}

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

calculate_d <- function(network_comm, names_col = "plant", abuns = NULL) {
  network_comm |>
    column_to_rownames(names_col) |>
    dfun(abuns = abuns)
}

extract_dstats <- function(dfun_rand_vals, d_stat_select) {
  dfun_rand_vals |>
    magrittr::extract(
      str_detect(names(dfun_rand_vals), d_stat_select)) |>
    purrr::list_c()
}

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

# Make the {bipartite} example dataset binary, then calculate a network
# index. This is to verify if some network indices always turn out 0.
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
