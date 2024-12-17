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

# Transpose the OTU community matrix
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
    bind_rows()
}

make_comm_for_network <- function(otu_clean, taxonomy_clean, plants) {
  otu_long_with_plants_and_taxonomy <-
    otu_clean |>
    pivot_longer(names_to = "sample", values_to = "abun", -otu_id) |>
    left_join(taxonomy_clean, by = "otu_id", relationship = "many-to-one") |>
    left_join(plants, by = "sample", relationship = "many-to-one")

}