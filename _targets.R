library(targets)
library(tarchetypes)
library(crew)

# パッケージと関数をロードする
source("R/packages.R")
source("R/functions.R")

# 並列処理の設定
# 重要：自分のパソコンのコア数を確認して、`workers`がそれを超えないようにすること
workers <- 16
if (parallel::detectCores() < workers) {
  stop("Number of workers exceeds number of cores")
}

tar_option_set(
  controller = crew_controller_local(workers = workers),
  workspace_on_error = TRUE
)

tar_plan(

  # データの読み込み ----

  # - ファイルの指定
  tar_file(
    otu_file,
    "data/otu_table.txt"
  ),
  tar_file(
    plants_file,
    "data/plants.csv"
  ),
  # - 植物データ
  plants = load_plants(plants_file),
  # - OTUデータ
  #   生データにはOTUのリード数と分類データが混合しているので、
  #   otu_taxonomyとotu_rawに分ける
  otu_raw_with_taxonomy = load_otu_raw(otu_file),
  otu_taxonomy = select(
    otu_raw_with_taxonomy, otu_id, taxonomy
  ),
  otu_raw = select(otu_raw_with_taxonomy, -taxonomy),

  # データクリーニング ----

  # 闘値を設定する（次のステップに使う）
  threshold_fraction = 0.001,
  # 低頻度OTUを除外する
  otu_clean = clean_otu_data(otu_raw, threshold_fraction),
  # 分類データを整える
  taxonomy_clean = clean_taxonomy_data(otu_taxonomy),
  # それぞれの分類階級がどれくらい同定できているのか確認する
  taxonomy_rank_count = count_tax_ranks(taxonomy_clean, otu_clean),

  # 希釈 ----

  # 準備：otuデータの列と行を入れ替える
  #  - otu_cleanではOTUが行、植物が列
  #  - otu_clean_transposedでは植物が行、OTUが列
  #  - グループ分けを設定する。各グループに６つの植物（行）が入っている
  tar_group_size(
    otu_clean_transposed,
    transpose_otu(otu_clean),
    6 # それぞれのグループに入っている植物の数
  ),
  # １グループずつiNEXTで希釈の解析を行う
  # 並列処理ができるので全部を一片に扱うより早い
  tar_target(
    i_next_res,
    i_next(otu_clean_transposed, nboot = 20),
    pattern = map(otu_clean_transposed),
    iteration = "list"
  ),
  asymptotic_richness = extract_asy_est(i_next_res),
  cor_res = calculate_diversity_corr(asymptotic_richness),

  # ネットワーク図 ----

  # ネットワーク図にデータを用意する
  # - 群集（community）データ。行が植物、列がOTUの分類群（共生菌）。#location_selectで地域選択
  network_comm_yoshimoto = make_comm_for_network(
    otu_clean, taxonomy_clean, plants, "otu_id", location_select = "吉本沖縄"
  ),
  
  network_comm_kumai = make_comm_for_network(
    otu_clean, taxonomy_clean, plants, "otu_id", location_select = "熊井沖縄"
  ),
  
  network_comm_zairai = make_comm_for_network(
    otu_clean, taxonomy_clean, plants, "otu_id", location_select = "在来種林"
  ),
  
  network_comm_akagi = make_comm_for_network(
    otu_clean, taxonomy_clean, plants, "otu_id", location_select = "アカギ純林"
  ),
  
  # - ネットワーク図のデータ
  network_graph_yoshimoto = make_network_graph(network_comm_yoshimoto),
  network_graph_kumai = make_network_graph(network_comm_kumai),
  network_graph_zairai = make_network_graph(network_comm_zairai),
  network_graph_akagi = make_network_graph(network_comm_akagi),
  
  # - 基本的な指標を計算する
  network_stats_obs_yoshimoto = calc_network_level(
    network_comm_yoshimoto,
    index = c("binary")
  ),
  network_stats_obs_kumai = calc_network_level(
    network_comm_kumai,
    index = c("binary")
  ),
  network_stats_obs_zairai = calc_network_level(
    network_comm_zairai,
    index = c("binary")
  ),
  network_stats_obs_akagi = calc_network_level(
    network_comm_akagi,
    index = c("binary")
  ),
  # 有意性検定の準備
  #  - curveballアルゴリズムで群集データを1万回いじる。
  #    図で元のでデータといじったデータを比較し、
  #    ランダムになるまで何回いじる必要があるのかを判断する。
  #    参照文献： https://docs.ropensci.org/canaper/articles/how-many-rand.html
  iter_sim_res_yoshimoto = cpr_iter_sim(
    comm = column_to_rownames(network_comm_yoshimoto, "plant"),
    null_model = "curveball", # 最も早いアルゴリズム
    n_iterations = 10000,
    thin = 10,
    seed = 123
  ),
  iter_sim_res_kumai = cpr_iter_sim(
    comm = column_to_rownames(network_comm_kumai, "plant"),
    null_model = "curveball", # 最も早いアルゴリズム
    n_iterations = 10000,
    thin = 10,
    seed = 123
  ),
  iter_sim_res_zairai = cpr_iter_sim(
    comm = column_to_rownames(network_comm_zairai, "plant"),
    null_model = "curveball", # 最も早いアルゴリズム
    n_iterations = 10000,
    thin = 10,
    seed = 123
  ),
  iter_sim_res_akagi = cpr_iter_sim(
    comm = column_to_rownames(network_comm_akagi, "plant"),
    null_model = "curveball", # 最も早いアルゴリズム
    n_iterations = 10000,
    thin = 10,
    seed = 123
  ),
  
  # - ランダムな群集を1000回シミュレーションする
  #   d 検定で効率よく並列処理ができるように100個のグループ（batch）に分ける
  #   各batchに10回シミュレーションする（合計1000回）
  tar_rep(
    random_comms_batched_yoshimoto,
    list(
      comm = randomize_single_comm(
        network_comm_yoshimoto,
        null_model = "curveball",
        n_iterations = 2000 # iter_sim_resの可視化によって2000回で十分であると判断
      )
    ),
    batches = 100,
    reps = 10,
    iteration = "list"
  ),
  tar_rep(
    random_comms_batched_kumai,
    list(
      comm = randomize_single_comm(
        network_comm_kumai,
        null_model = "curveball",
        n_iterations = 2000 # iter_sim_resの可視化によって2000回で十分であると判断
      )
    ),
    batches = 100,
    reps = 10,
    iteration = "list"
  ),
  tar_rep(
    random_comms_batched_zairai,
    list(
      comm = randomize_single_comm(
        network_comm_zairai,
        null_model = "curveball",
        n_iterations = 2000 # iter_sim_resの可視化によって2000回で十分であると判断
      )
    ),
    batches = 100,
    reps = 10,
    iteration = "list"
  ),
  tar_rep(
    random_comms_batched_akagi,
    list(
      comm = randomize_single_comm(
        network_comm_akagi,
        null_model = "curveball",
        n_iterations = 2000 # iter_sim_resの可視化によって2000回で十分であると判断
      )
    ),
    batches = 100,
    reps = 10,
    iteration = "list"
  ),
  # - 検定する指標を設定する
  #   weighted NODFとH2はバイナリーデータで必ずゼロになるので、ここでは使わない
  network_indices = c("NODF", "connectance"),
  # - 有意性検定を行う（ネットワーク全体の指標）
  tar_target(
    random_comms_list_yoshimoto,
    purrr::flatten(random_comms_batched_yoshimoto) |> map(1)
  ),
  tar_target(
    network_stats_yoshimoto,
    run_rand_test(
      network_comm_yoshimoto,
      random_comms_list_yoshimoto,
      index = network_indices
    ),
    pattern = map(network_indices)
  ),
  tar_target(
    random_comms_list_kumai,
    purrr::flatten(random_comms_batched_kumai) |> map(1)
  ),
  tar_target(
    network_stats_kumai,
    run_rand_test(
      network_comm_kumai,
      random_comms_list_kumai,
      index = network_indices
    ),
    pattern = map(network_indices)
  ),
  tar_target(
    random_comms_list_zairai,
    purrr::flatten(random_comms_batched_zairai) |> map(1)
  ),
  tar_target(
    network_stats_zairai,
    run_rand_test(
      network_comm_zairai,
      random_comms_list_zairai,
      index = network_indices
    ),
    pattern = map(network_indices)
  ),
  tar_target(
    random_comms_list_akagi,
    purrr::flatten(random_comms_batched_akagi) |> map(1)
  ),
  tar_target(
    network_stats_akagi,
    run_rand_test(
      network_comm_akagi,
      random_comms_list_akagi,
      index = network_indices
    ),
    pattern = map(network_indices)
  ),
  
  # - 種特殊性の有意性検定を行う（種ごとの指標）
  dfun_obs_res_yoshimoto = calculate_d(network_comm_yoshimoto),
  tar_rep2(
    dfun_rand_vals_yoshimoto,
    dfun(random_comms_batched_yoshimoto$comm),
    random_comms_batched_yoshimoto
  ),
  plant_names_yoshimoto = pull(network_comm_yoshimoto, "plant"),
  d_indices = c("dprime", "d", "dmin", "dmax"),
  tar_target(
    d_stats_yoshimoto,
    calculate_p_val_dstats(
      dfun_obs_res = dfun_obs_res_yoshimoto,
      dfun_rand_vals = dfun_rand_vals_yoshimoto,
      d_stat_select = d_indices,
      species_select = plant_names_yoshimoto
    ),
    pattern = cross(d_indices, plant_names_yoshimoto)
  ),
  
  
  dfun_obs_res_kumai = calculate_d(network_comm_kumai),
  tar_rep2(
    dfun_rand_vals_kumai,
    dfun(random_comms_batched_kumai$comm),
    random_comms_batched_kumai
  ),
  plant_names_kumai = pull(network_comm_kumai, "plant"),
  tar_target(
    d_stats_kumai,
    calculate_p_val_dstats(
      dfun_obs_res = dfun_obs_res_kumai,
      dfun_rand_vals = dfun_rand_vals_kumai,
      d_stat_select = d_indices,
      species_select = plant_names_kumai
    ),
    pattern = cross(d_indices, plant_names_kumai)
  ),
  
  
  dfun_obs_res_zairai = calculate_d(network_comm_zairai),
  tar_rep2(
    dfun_rand_vals_zairai,
    dfun(random_comms_batched_zairai$comm),
    random_comms_batched_zairai
  ),
  plant_names_zairai = pull(network_comm_zairai, "plant"),
  tar_target(
    d_stats_zairai,
    calculate_p_val_dstats(
      dfun_obs_res = dfun_obs_res_zairai,
      dfun_rand_vals = dfun_rand_vals_zairai,
      d_stat_select = d_indices,
      species_select = plant_names_zairai
    ),
    pattern = cross(d_indices, plant_names_zairai)
  ),
  
  
  dfun_obs_res_akagi = calculate_d(network_comm_akagi),
  tar_rep2(
    dfun_rand_vals_akagi,
    dfun(random_comms_batched_akagi$comm),
    random_comms_batched_akagi
  ),
  plant_names_akagi = pull(network_comm_akagi, "plant"),
  tar_target(
    d_stats_akagi,
    calculate_p_val_dstats(
      dfun_obs_res = dfun_obs_res_akagi,
      dfun_rand_vals = dfun_rand_vals_akagi,
      d_stat_select = d_indices,
      species_select = plant_names_akagi
    ),
    pattern = cross(d_indices, plant_names_akagi)
  ),
  
  #真菌から植物の種特殊性
  #行列の入れ替え
  Tnetwork_comm_yoshimoto <- t(network_comm_yoshimoto),
  colnames(Tnetwork_comm_yoshimoto) <- as.character(Tnetwork_comm_yoshimoto[1,]),
  Tnetwork_comm_yoshimoto <- Tnetwork_comm_yoshimoto[-1,],
   
  Tnetwork_comm_kumai <- t(network_comm_kumai),
  colnames(Tnetwork_comm_kumai) <- as.character(Tnetwork_comm_kumai[1,]),
  Tnetwork_comm_kumai <- Tnetwork_comm_kumai[-1,],
  
  Tnetwork_comm_zairai <- t(network_comm_zairai),
  colnames(Tnetwork_comm_zairai) <- as.character(Tnetwork_comm_zairai[1,]),
  Tnetwork_comm_zairai <- Tnetwork_comm_zairai[-1,],
  
  Tnetwork_comm_akagi <- t(network_comm_akagi),
  colnames(Tnetwork_comm_akagi) <- as.character(Tnetwork_comm_akagi[1,]),
  Tnetwork_comm_akagi <- Tnetwork_comm_akagi[-1,],
  

  # - リポートを書く
  tar_quarto(
    network_report_yoshimoto,
    path = "report.qmd",
    quiet = FALSE
  ),

 tar_quarto(
  network_report_kumai,
  path = "report.qmd",
  quiet = FALSE
),


 tar_quarto(
  network_report_zairai,
  path = "report.qmd",
  quiet = FALSE
),

 tar_quarto(
  network_report_akagi,
  path = "report.qmd",
  quiet = FALSE
)
)
