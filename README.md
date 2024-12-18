
# seitai-yomitoku-targets

このリポジトリは、「DNA情報で生態系を読み解く」（東樹 宏和 著）で紹介されている解析を再現するためのコードを提供します。

This repository provides code for replicating (some) of the analyses in "DNA joho de seitai wo yomitoku" by Hirokazu Toju.

## コードの実行

このプロジェクトは[targets](https://github.com/ropensci/targets)でワークフロー管理を行なっています。

下の準備が終わったら、Rで以下のコマンドを実行するとワークフローが自動的に実行されます：

```
targets::tar_make()
```

## How to open in RStudio | RStudioでの開き方

1. RStudio の File メニューから New Project を選択して、次に Version Control を選択します。  
![](docs/r-project-wizard-1.png)
2. Git を選択します。  
![](docs/r-project-wizard-2.png)
3. Repository URLに <https://github.com/joelnitta/seitai-yomitoku/>、Project directory name に`seitai-yomitoku`、Create project as subdirectory of に任意のフォルダー（Desktopがおすすめ）を入力します。  
![](docs/r-project-wizard-3.png)
4. Create project ボタンを押します。
5. データをダウンロードする

## Data | データ

使用するデータファイルは以下の2つです：

There are two data files:

- `otu_table.txt`  
- `plants.csv`

これらのファイルは[Google Drive](https://drive.google.com/drive/folders/1obxIYrq2isURX79c0Skm5TGeO5KhbPrq?usp=drive_link)（現在アクセスが限定されています）からダウンロードしてください。  
ダウンロード後、`data/`ディレクトリ内に配置してください。

Download these from [Google Drive](https://drive.google.com/drive/folders/1obxIYrq2isURX79c0Skm5TGeO5KhbPrq?usp=drive_link) (currently, access is restricted) and put them in the `data/` folder.


## ライセンス | License

seitai-yomitoku is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

seitai-yomitokuは、[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)の下でライセンスされています。