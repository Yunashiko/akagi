library(targets)
library(tarchetypes)
library(crew)

# パッケージと関数をロードする
source("R/packages.R")
source("R/functions.R")

# 並列処理の設定
# 重要：自分のパソコンのコア数を確認して、`workers`がそれを超えないようにすること
workers <- 6
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
  taxonomy_rank_count = count_tax_ranks(taxonomy_clean),

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
  # - 群集（community）データ。行が植物、列がOTUの分類群（共生菌）。
  network_comm = make_comm_for_network(
    otu_clean, taxonomy_clean, plants, "family"
  ),
  # - ネットワーク図のデータ
  network_graph = make_network_graph(network_comm),
  # - 基本的な指標を計算する
  network_stats_obs = calc_network_level(
    network_comm,
    index = c("binary")
  ),
  # 有意性検定の準備
  #  - curveballアルゴリズムで群集データを1万回いじる。
  #    図で元のでデータといじったデータを比較し、
  #    ランダムになるまで何回いじる必要があるのかを判断する。
  #    参照文献： https://docs.ropensci.org/canaper/articles/how-many-rand.html
  iter_sim_res = cpr_iter_sim(
    comm = column_to_rownames(network_comm, "plant"),
    null_model = "curveball", # 最も早いアルゴリズム
    n_iterations = 10000,
    thin = 10,
    seed = 123
  ),
  # - ランダムな群集を1000回シミュレーションする
  #   d 検定で効率よく並列処理ができるように100個のグループ（batch）に分ける
  #   各batchに10回シミュレーションする（合計1000回）
  tar_rep(
    random_comms_batched,
    list(
      comm = randomize_single_comm(
        network_comm,
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
  network_indices = c("NODF"),
  # - 有意性検定を行う（ネットワーク全体の指標）
  tar_target(
    random_comms_list,
    purrr::flatten(random_comms_batched) |> map(1)
  ),
  tar_target(
    network_stats,
    run_rand_test(
      network_comm,
      random_comms_list,
      index = network_indices
    ),
    pattern = map(network_indices)
  ),
  # - 種特殊性の有意性検定を行う（種ごとの指標）
  dfun_obs_res = calculate_d(network_comm),
  tar_rep2(
    dfun_rand_vals,
    dfun(random_comms_batched$comm),
    random_comms_batched
  ),
  plant_names = pull(network_comm, "plant"),
  d_indices = c("dprime", "d", "dmin", "dmax"),
  tar_target(
    d_stats,
    calculate_p_val_dstats(
      dfun_obs_res = dfun_obs_res,
      dfun_rand_vals = dfun_rand_vals,
      d_stat_select = d_indices,
      species_select = plant_names
    ),
    pattern = cross(d_indices, plant_names)
  ),
  # - リポートを書く
  tar_quarto(
    network_report,
    path = "report.qmd",
    quiet = FALSE
  )
)
