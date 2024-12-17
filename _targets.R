library(targets)
library(tarchetypes)
library(crew)

# パッケージと関数をロードする
source("R/packages.R")
source("R/functions.R")

# 並列処理の設定
# 重要：自分のパソコンのコア数を確認して、`workers`がそれを超えないようにすること
tar_option_set(
  controller = crew_controller_local(workers = 6)
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
  # 生データにはOTUのリード数と分類データが混合しているので、分ける
  otu_raw_with_taxonomy = load_otu_raw(otu_file),
  otu_taxonomy = select(
    otu_raw_with_taxonomy, otu_id, taxonomy),
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
  cor_res = check_diversity_correlation(asymptotic_richness)
)
