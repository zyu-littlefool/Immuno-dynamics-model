# ============================================================
# BCR 14/21/28 天联合谱系轨迹图
#
# 功能：
# 1) 自动读取同一文件夹中的 CSV/TSV/TXT。
# 2) 从文件名识别实验组以及 14、21、28 天。
#    同时支持 14/21/28、模拟步数 1400/2100/2800，
#    以及 MATLAB 下标形式 1401/2101/2801。
# 3) 每个时间点先合并完全相同的 char，再选取前 top_n_per_time 个 clone。
# 4) 使用 char 的 Hamming 距离，把后一个时间点的 clone 连接到
#    最相似的较早 clone；一个父 clone 可以产生多个后代分支。
# 5) 将 14、21、28 天放在三个同心时间环中，每个实验组单独输出一张图。
# 6) 输出节点表和父子匹配表，便于检查“合并谱系”的依据。
#
# 重要说明：
# 仅凭 char 相似性得到的是“推断谱系轨迹”，不是严格的 BCR 系谱证明。
# 如果有 V/J 基因、CDR3、克隆型编号或 germline 序列，应当先按这些信息
# 划分 clonotype，再在 clonotype 内进行 Hamming 距离匹配。
# ============================================================

# -----------------------------
# 0. 包检查与加载
# -----------------------------
cran_packages <- c(
  "dplyr", "stringdist", "ggplot2", "readr", "scales"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("正在安装缺失的包：", pkg)
    install.packages(pkg)
  } else {
    message("已安装：", pkg)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringdist)
  library(ggplot2)
  library(readr)
  library(scales)
})

# -----------------------------
# 1. 文件读取与文件名解析
# -----------------------------

infer_delimiter <- function(file_path, encoding = "UTF-8") {
  header <- readLines(
    file_path,
    n = 1,
    warn = FALSE,
    encoding = encoding
  )
  
  if (length(header) == 0) {
    stop("文件为空：", file_path)
  }
  
  if (grepl(",", header, fixed = TRUE)) {
    return(",")
  }
  if (grepl("\t", header, fixed = TRUE)) {
    return("\t")
  }
  if (grepl(";", header, fixed = TRUE)) {
    return(";")
  }
  
  stop("无法判断文件分隔符：", file_path)
}

read_bcr_csv <- function(file_path, encoding = "UTF-8") {
  delim <- infer_delimiter(file_path, encoding = encoding)
  
  df <- readr::read_delim(
    file = file_path,
    delim = delim,
    locale = readr::locale(encoding = encoding),
    col_types = readr::cols(
      PatternID    = readr::col_character(),
      Count        = readr::col_double(),
      Frequency    = readr::col_double(),
      E_change_wt  = readr::col_double(),
      E_change_var = readr::col_double(),
      char         = readr::col_character(),
      .default     = readr::col_guess()
    ),
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  required_cols <- c(
    "PatternID", "Count", "E_change_var", "char"
  )
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "文件 ", basename(file_path),
      " 缺少必要列：",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      PatternID = as.character(.data$PatternID),
      Count = as.numeric(.data$Count),
      E_change_var = as.numeric(.data$E_change_var),
      char = gsub("\\s+", "", as.character(.data$char))
    )
}

# 为某一天生成可识别的文件名正则表达式。
# 示例：
# fast_gc_14.csv
# fast_gc_day14.csv
# fast_gc_14d.csv
# fast_gc_1400.csv
# Gene_fast_at_1401_GC_sub_relmod1_seed11_mIter34.csv
day_filename_patterns <- function(day_value) {
  simulation_step <- day_value * 100
  # MATLAB 数组下标通常比模拟步数多 1：
  # 第 14/21/28 天也可能写成 1401/2101/2801。
  simulation_index <- simulation_step + 1
  
  c(
    # 优先匹配并整体删除 at_1401、at_2101、at_2801，
    # 这样解析后的 Group 中不会残留单独的 "at"。
    paste0(
      "(^|[^[:alnum:]])at[_ .-]*(?:",
      simulation_step,
      "|",
      simulation_index,
      ")([^0-9]|$)"
    ),
    paste0(
      "(^|[^0-9])(?:",
      simulation_step,
      "|",
      simulation_index,
      ")",
      "([^0-9]|$)"
    ),
    paste0(
      "(^|[^[:alnum:]])(?:day|d)[_ .-]*",
      day_value,
      "([^0-9]|$)"
    ),
    paste0(
      "(^|[^0-9])",
      day_value,
      "[_ .-]*(?:day|d|天)([^[:alnum:]]|$)"
    ),
    paste0(
      "(^|[^0-9])",
      day_value,
      "([^0-9]|$)"
    )
  )
}

parse_group_and_day <- function(file_path,
                                days = c(14, 21, 28)) {
  stem <- tools::file_path_sans_ext(basename(file_path))
  stem_lower <- tolower(stem)
  
  hit_days <- numeric(0)
  hit_patterns <- character(0)
  
  for (day_value in days) {
    patterns <- day_filename_patterns(day_value)
    
    for (pattern in patterns) {
      match_position <- regexpr(
        pattern,
        stem_lower,
        perl = TRUE
      )
      
      if (match_position[1] > 0) {
        hit_days <- c(hit_days, day_value)
        hit_patterns <- c(hit_patterns, pattern)
        break
      }
    }
  }
  
  unique_days <- unique(hit_days)
  
  if (length(unique_days) == 0) {
    return(data.frame(
      Sample = stem,
      Group = NA_character_,
      Day = NA_real_,
      ParseStatus = paste0(
        "未识别到 14/21/28、1400/2100/2800",
        " 或 1401/2101/2801"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(unique_days) > 1) {
    return(data.frame(
      Sample = stem,
      Group = NA_character_,
      Day = NA_real_,
      ParseStatus = paste0(
        "文件名同时包含多个时间点：",
        paste(unique_days, collapse = ",")
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  day_value <- unique_days[1]
  pattern_to_remove <- hit_patterns[
    which(hit_days == day_value)[1]
  ]
  
  group_name <- gsub(
    pattern_to_remove,
    "_",
    stem,
    perl = TRUE,
    ignore.case = TRUE
  )
  group_name <- gsub("[[:space:].-]+", "_", group_name)
  group_name <- gsub("_+", "_", group_name)
  group_name <- gsub("^_+|_+$", "", group_name)
  
  if (group_name == "") {
    group_name <- "group"
  }
  
  data.frame(
    Sample = stem,
    Group = group_name,
    Day = day_value,
    ParseStatus = "OK",
    stringsAsFactors = FALSE
  )
}

build_file_manifest <- function(csv_dir,
                                days = c(14, 21, 28),
                                pattern = "\\.(csv|tsv|txt)$",
                                strict_filename_parsing = FALSE) {
  files <- list.files(
    path = csv_dir,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = FALSE
  )
  
  if (length(files) == 0) {
    stop(
      "文件夹中没有找到 CSV/TSV/TXT：",
      normalizePath(csv_dir, mustWork = FALSE)
    )
  }
  
  parsed_rows <- lapply(
    files,
    parse_group_and_day,
    days = days
  )
  
  manifest <- bind_rows(parsed_rows)
  manifest$File <- files
  
  bad_rows <- manifest$ParseStatus != "OK"
  
  if (any(bad_rows)) {
    bad_text <- paste0(
      basename(manifest$File[bad_rows]),
      " [",
      manifest$ParseStatus[bad_rows],
      "]"
    )
    
    if (strict_filename_parsing) {
      stop(
        "以下文件名无法解析：\n",
        paste(bad_text, collapse = "\n")
      )
    } else {
      warning(
        "以下文件不会参加绘图：\n",
        paste(bad_text, collapse = "\n")
      )
      manifest <- manifest[!bad_rows, , drop = FALSE]
    }
  }
  
  if (nrow(manifest) == 0) {
    stop("没有成功识别任何 14/21/28 天数据文件。")
  }
  
  manifest %>%
    arrange(.data$Group, .data$Day, .data$Sample) %>%
    select(
      .data$File,
      .data$Sample,
      .data$Group,
      .data$Day,
      .data$ParseStatus
    )
}

read_manifest_data <- function(manifest,
                               encoding = "UTF-8") {
  data_rows <- vector("list", nrow(manifest))
  
  for (row_id in seq_len(nrow(manifest))) {
    one_df <- read_bcr_csv(
      manifest$File[row_id],
      encoding = encoding
    )
    
    one_df$SourceFile <- basename(
      manifest$File[row_id]
    )
    one_df$Sample <- manifest$Sample[row_id]
    one_df$Group <- manifest$Group[row_id]
    one_df$Day <- manifest$Day[row_id]
    
    data_rows[[row_id]] <- one_df
  }
  
  bind_rows(data_rows)
}

# -----------------------------
# 2. 每个时间点筛选与重复序列合并
# -----------------------------

weighted_mean_safe <- function(value, weight) {
  value <- as.numeric(value)
  weight <- as.numeric(weight)
  
  keep <- is.finite(value)
  value <- value[keep]
  weight <- weight[keep]
  
  if (length(value) == 0) {
    return(NA_real_)
  }
  
  weight[!is.finite(weight) | weight < 0] <- 0
  
  if (sum(weight) <= 0) {
    return(mean(value))
  }
  
  sum(value * weight) / sum(weight)
}

classify_count_for_plot <- function(count) {
  count <- as.numeric(count)
  
  dplyr::case_when(
    is.na(count) | count <= 1 ~ "1",
    count == 2 ~ "2",
    count == 3 ~ "3",
    count == 4 ~ "4",
    count >= 5 ~ ">=5",
    TRUE ~ "1"
  )
}

prepare_group_nodes <- function(group_df,
                                group_name,
                                days = c(14, 21, 28),
                                top_n_per_time = 40) {
  group_df <- group_df %>%
    filter(
      .data$Day %in% days,
      !is.na(.data$char),
      .data$char != ""
    ) %>%
    mutate(
      Count = as.numeric(.data$Count),
      E_change_var = as.numeric(.data$E_change_var)
    )
  
  available_days <- sort(unique(group_df$Day))
  missing_days <- setdiff(days, available_days)
  
  if (length(missing_days) > 0) {
    stop(
      "组 ", group_name,
      " 缺少时间点：",
      paste(missing_days, collapse = ", ")
    )
  }
  
  # 同一个时间点中完全相同的 char 先合并。
  # Count 求和；E_change_var 按 Count 加权平均。
  collapsed <- group_df %>%
    group_by(.data$Day, .data$char) %>%
    summarise(
      PatternIDs = paste(
        unique(.data$PatternID),
        collapse = ";"
      ),
      E_change_var = weighted_mean_safe(
        .data$E_change_var,
        .data$Count
      ),
      Count = sum(.data$Count, na.rm = TRUE),
      SourceFiles = paste(
        unique(.data$SourceFile),
        collapse = ";"
      ),
      .groups = "drop"
    )
  
  selected <- collapsed %>%
    group_by(.data$Day) %>%
    arrange(
      desc(.data$E_change_var),
      desc(.data$Count),
      .data$char,
      .by_group = TRUE
    ) %>%
    slice_head(n = top_n_per_time) %>%
    mutate(
      RankWithinDay = row_number(),
      NodeID = sprintf(
        "D%02d_N%03d",
        .data$Day,
        .data$RankWithinDay
      )
    ) %>%
    ungroup() %>%
    mutate(
      Group = group_name,
      CountGroup = factor(
        classify_count_for_plot(.data$Count),
        levels = c("1", "2", "3", "4", ">=5")
      )
    )
  
  sequence_lengths <- nchar(selected$char)
  
  if (length(unique(sequence_lengths)) != 1) {
    length_table <- table(sequence_lengths)
    
    stop(
      "组 ", group_name,
      " 的 char 长度不一致，不能计算 Hamming 距离。长度分布：",
      paste(
        names(length_table),
        length_table,
        sep = "=",
        collapse = "; "
      )
    )
  }
  
  selected %>%
    arrange(
      match(.data$Day, days),
      .data$RankWithinDay
    )
}

# -----------------------------
# 3. 跨时间点推断父子谱系
# -----------------------------

calculate_hamming_threshold <- function(sequence_length,
                                        max_hamming_distance = NULL,
                                        max_hamming_fraction = 0.10) {
  if (!is.null(max_hamming_distance)) {
    max_hamming_distance <- as.numeric(
      max_hamming_distance
    )
    
    if (length(max_hamming_distance) != 1 ||
        !is.finite(max_hamming_distance) ||
        max_hamming_distance < 0) {
      stop(
        "max_hamming_distance 必须是 NULL 或非负有限数值。"
      )
    }
    
    return(max_hamming_distance)
  }
  
  max_hamming_fraction <- as.numeric(
    max_hamming_fraction
  )
  
  if (length(max_hamming_fraction) != 1 ||
      !is.finite(max_hamming_fraction) ||
      max_hamming_fraction < 0 ||
      max_hamming_fraction > 1) {
    stop(
      "max_hamming_fraction 必须位于 0 到 1 之间。"
    )
  }
  
  max(
    1,
    floor(sequence_length * max_hamming_fraction)
  )
}

infer_temporal_edges <- function(nodes,
                                 days = c(14, 21, 28),
                                 max_hamming_distance = NULL,
                                 max_hamming_fraction = 0.10,
                                 root_id = "ROOT") {
  sequence_length <- unique(nchar(nodes$char))
  
  if (length(sequence_length) != 1) {
    stop("infer_temporal_edges 要求所有 char 长度一致。")
  }
  
  threshold <- calculate_hamming_threshold(
    sequence_length = sequence_length,
    max_hamming_distance = max_hamming_distance,
    max_hamming_fraction = max_hamming_fraction
  )
  
  root_char <- paste0(
    rep("0", sequence_length),
    collapse = ""
  )
  
  edge_rows <- list()
  edge_counter <- 0
  
  for (day_value in days) {
    children <- nodes[
      nodes$Day == day_value,
      ,
      drop = FALSE
    ]
    
    for (child_row in seq_len(nrow(children))) {
      child <- children[
        child_row,
        ,
        drop = FALSE
      ]
      
      parent_id <- root_id
      parent_day <- 0
      parent_distance <- stringdist::stringdist(
        root_char,
        child$char,
        method = "hamming"
      )
      match_type <- if (day_value == min(days)) {
        "root_start"
      } else {
        "new_lineage"
      }
      
      earlier_days <- rev(
        days[days < day_value]
      )
      
      # 优先在紧邻的上一个时间点寻找父节点；
      # 找不到合格节点时，再回退到更早的时间点。
      if (length(earlier_days) > 0) {
        for (candidate_day in earlier_days) {
          candidates <- nodes[
            nodes$Day == candidate_day,
            ,
            drop = FALSE
          ]
          
          if (nrow(candidates) == 0) {
            next
          }
          
          candidate_distance <- as.numeric(
            stringdist::stringdist(
              child$char,
              candidates$char,
              method = "hamming"
            )
          )
          
          eligible <- which(
            candidate_distance <= threshold
          )
          
          if (length(eligible) == 0) {
            next
          }
          
          candidate_rank <- data.frame(
            RowID = eligible,
            Distance = candidate_distance[eligible],
            CountSort = -ifelse(
              is.finite(candidates$Count[eligible]),
              candidates$Count[eligible],
              0
            ),
            EvarSort = -ifelse(
              is.finite(
                candidates$E_change_var[eligible]
              ),
              candidates$E_change_var[eligible],
              -Inf
            ),
            NodeID = candidates$NodeID[eligible],
            stringsAsFactors = FALSE
          )
          
          candidate_rank <- candidate_rank[
            order(
              candidate_rank$Distance,
              candidate_rank$CountSort,
              candidate_rank$EvarSort,
              candidate_rank$NodeID
            ),
            ,
            drop = FALSE
          ]
          
          winner_row <- candidate_rank$RowID[1]
          parent_id <- candidates$NodeID[winner_row]
          parent_day <- candidates$Day[winner_row]
          parent_distance <- candidate_distance[
            winner_row
          ]
          match_type <- if (
            candidate_day ==
            max(days[days < day_value])
          ) {
            "matched_previous_time"
          } else {
            "matched_earlier_time"
          }
          
          break
        }
      }
      
      edge_counter <- edge_counter + 1
      edge_rows[[edge_counter]] <- data.frame(
        EdgeID = sprintf(
          "E%04d",
          edge_counter
        ),
        ParentID = parent_id,
        ChildID = child$NodeID,
        ParentDay = parent_day,
        ChildDay = child$Day,
        HammingDistance = parent_distance,
        Similarity = 1 -
          parent_distance / sequence_length,
        MatchThreshold = threshold,
        MatchType = match_type,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(edge_rows)
}

# -----------------------------
# 4. 同心时间环布局
# -----------------------------

order_children_by_sequence <- function(child_ids,
                                       plot_nodes) {
  child_ids <- unique(as.character(child_ids))
  
  if (length(child_ids) <= 1) {
    return(child_ids)
  }
  
  row_index <- match(
    child_ids,
    plot_nodes$NodeID
  )
  child_char <- plot_nodes$char[row_index]
  
  if (any(is.na(child_char)) ||
      length(unique(nchar(child_char))) != 1) {
    return(sort(child_ids))
  }
  
  if (length(child_ids) == 2) {
    return(
      child_ids[
        order(child_char, child_ids)
      ]
    )
  }
  
  distance_matrix <- as.matrix(
    stringdist::stringdistmatrix(
      child_char,
      child_char,
      method = "hamming"
    )
  )
  
  if (all(distance_matrix == 0)) {
    return(sort(child_ids))
  }
  
  cluster_result <- hclust(
    as.dist(distance_matrix),
    method = "average"
  )
  
  child_ids[cluster_result$order]
}

make_radial_time_layout <- function(nodes,
                                    edges,
                                    days = c(14, 21, 28),
                                    root_id = "ROOT",
                                    start_angle = pi / 2,
                                    top_label_gap = 0) {
  root_row <- data.frame(
    Group = unique(nodes$Group)[1],
    Day = 0,
    char = paste0(
      rep("0", unique(nchar(nodes$char))),
      collapse = ""
    ),
    PatternIDs = root_id,
    Count = NA_real_,
    E_change_var = 0,
    SourceFiles = NA_character_,
    RankWithinDay = 0,
    NodeID = root_id,
    CountGroup = factor(
      NA_character_,
      levels = c("1", "2", "3", "4", ">=5")
    ),
    stringsAsFactors = FALSE
  )
  
  plot_nodes <- bind_rows(
    root_row,
    nodes
  )
  
  children_by_parent <- split(
    edges$ChildID,
    edges$ParentID
  )
  
  leaf_cache <- new.env(
    hash = TRUE,
    parent = emptyenv()
  )
  
  get_descendant_leaves <- function(node_id) {
    cache_name <- paste0("N_", node_id)
    
    if (exists(
      cache_name,
      envir = leaf_cache,
      inherits = FALSE
    )) {
      return(
        get(
          cache_name,
          envir = leaf_cache,
          inherits = FALSE
        )
      )
    }
    
    children <- children_by_parent[[node_id]]
    
    if (is.null(children) ||
        length(children) == 0) {
      answer <- node_id
    } else {
      ordered_children <- order_children_by_sequence(
        children,
        plot_nodes
      )
      answer <- unlist(
        lapply(
          ordered_children,
          get_descendant_leaves
        ),
        use.names = FALSE
      )
    }
    
    assign(
      cache_name,
      answer,
      envir = leaf_cache
    )
    answer
  }
  
  leaf_order <- unique(
    get_descendant_leaves(root_id)
  )
  
  if (length(leaf_order) == 0) {
    stop("没有可用于布局的末端节点。")
  }
  
  # 在圆形图正上方保留一个小缺口，防止时间标签遮住节点。
  leaf_angle <- start_angle + top_label_gap / 2 +
    seq(
      0,
      2 * pi - top_label_gap,
      length.out = length(leaf_order) + 1
    )[
      seq_along(leaf_order)
    ]
  
  names(leaf_angle) <- leaf_order
  
  node_angle <- vapply(
    plot_nodes$NodeID,
    function(node_id) {
      descendants <- get_descendant_leaves(
        node_id
      )
      mean(leaf_angle[descendants])
    },
    numeric(1)
  )
  
  time_radii <- setNames(
    c(1.4, 2.2, 3.0),
    as.character(days)
  )
  
  node_radius <- ifelse(
    plot_nodes$Day == 0,
    0,
    unname(
      time_radii[
        as.character(plot_nodes$Day)
      ]
    )
  )
  
  plot_nodes$Angle <- node_angle
  plot_nodes$Radius <- node_radius
  plot_nodes$X <- node_radius * cos(node_angle)
  plot_nodes$Y <- node_radius * sin(node_angle)
  
  has_children <- plot_nodes$NodeID %in%
    edges$ParentID
  
  plot_nodes$LineageStatus <- case_when(
    plot_nodes$NodeID == root_id ~ "Root",
    plot_nodes$Day < max(days) &
      !has_children ~ "Ended before final time",
    has_children ~ "Continuing",
    TRUE ~ "Observed at final time"
  )
  
  plot_nodes
}

make_edge_path_data <- function(edges,
                                layout_nodes,
                                points_per_edge = 30,
                                root_id = "ROOT") {
  path_rows <- vector(
    "list",
    nrow(edges)
  )
  
  for (edge_row in seq_len(nrow(edges))) {
    edge <- edges[
      edge_row,
      ,
      drop = FALSE
    ]
    
    parent <- layout_nodes[
      layout_nodes$NodeID == edge$ParentID,
      ,
      drop = FALSE
    ]
    child <- layout_nodes[
      layout_nodes$NodeID == edge$ChildID,
      ,
      drop = FALSE
    ]
    
    if (nrow(parent) != 1 ||
        nrow(child) != 1) {
      stop(
        "无法为边 ",
        edge$EdgeID,
        " 找到唯一的父节点和子节点。"
      )
    }
    
    t_value <- seq(
      0,
      1,
      length.out = points_per_edge
    )
    
    parent_angle <- parent$Angle
    child_angle <- child$Angle
    
    if (parent$NodeID == root_id) {
      parent_angle <- child_angle
    }
    
    # smoothstep 让分支从父节点自然弯向子节点。
    smooth_value <- t_value^2 *
      (3 - 2 * t_value)
    
    path_radius <- parent$Radius +
      (child$Radius - parent$Radius) *
      t_value
    
    path_angle <- parent_angle +
      (child_angle - parent_angle) *
      smooth_value
    
    path_rows[[edge_row]] <- data.frame(
      EdgeID = edge$EdgeID,
      PathOrder = seq_along(t_value),
      X = path_radius * cos(path_angle),
      Y = path_radius * sin(path_angle),
      MatchType = edge$MatchType,
      HammingDistance =
        edge$HammingDistance,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(path_rows)
}

make_time_ring_data <- function(
    days = c(14, 21, 28),
    points_per_ring = 361) {
  
  # 必须与节点使用的时间半径完全一致
  time_radii <- c(1.4, 2.2, 3.0)
  
  ring_rows <- vector(
    "list",
    length(days)
  )
  
  for (ring_id in seq_along(days)) {
    theta <- seq(
      0,
      2 * pi,
      length.out = points_per_ring
    )
    
    ring_radius <- time_radii[ring_id]
    
    ring_rows[[ring_id]] <- data.frame(
      Day = days[ring_id],
      RingID = ring_id,
      X = ring_radius * cos(theta),
      Y = ring_radius * sin(theta)
    )
  }
  
  bind_rows(ring_rows)
}

# -----------------------------
# 5. 绘图
# -----------------------------

get_global_evar_limits <- function(node_list) {
  values <- unlist(
    lapply(
      node_list,
      function(x) x$E_change_var
    )
  )
  values <- values[is.finite(values)]
  
  if (length(values) == 0) {
    return(c(-1, 1))
  }
  
  lower <- min(c(values, 0))
  upper <- max(c(values, 1))
  
  if (lower == upper) {
    lower <- lower - 0.1
    upper <- upper + 0.1
  }
  
  c(lower, upper)
}

plot_group_temporal_tree <- function(group_name,
                                     layout_nodes,
                                     edges,
                                     days = c(14, 21, 28),
                                     evar_limits = c(-1, 1),
                                     show_node_labels = FALSE) {
  edge_paths <- make_edge_path_data(
    edges = edges,
    layout_nodes = layout_nodes
  )
  
  ring_data <- make_time_ring_data(days)
  
  clone_nodes <- layout_nodes %>%
    filter(.data$Day > 0)
  
  root_node <- layout_nodes %>%
    filter(.data$Day == 0)
  
  ended_nodes <- clone_nodes %>%
    filter(
      .data$LineageStatus ==
        "Ended before final time"
    )
  
 
  count_levels <- c(
    "1", "2", "3", "4", ">=5"
  )
  count_labels <- c(
    "1", "2", "3", "4", "≥5"
  )
  count_sizes <- c(
    "1" = 2.8,
    "2" = 3.8,
    "3" = 5.0,
    "4" = 6.2,
    ">=5" = 7.4
  )
  
  max_radius <- length(days)
  plotting_limit <- max_radius + 0.65
  
  p <- ggplot() +
    geom_path(
      data = ring_data,
      aes(
        x = .data$X,
        y = .data$Y,
        group = .data$RingID
      ),
      color = "#53656D",
      linewidth = 0.55,
      linetype = "dotted"
    ) +
    geom_path(
      data = edge_paths,
      aes(
        x = .data$X,
        y = .data$Y,
        group = .data$EdgeID
      ),
      color = "#53656D",
      linewidth = 0.75,
      lineend = "round"
    ) +
    geom_point(
      data = clone_nodes,
      aes(
        x = .data$X,
        y = .data$Y,
        size = .data$CountGroup,
        color = .data$E_change_var
      ),
      alpha = 0.96,
      stroke = 0.5
    ) +
    geom_point(
      data = ended_nodes,
      aes(
        x = .data$X,
        y = .data$Y,
        size = .data$CountGroup
      ),
      shape = 1,
      color = "#657983",
      stroke = 1.1,
      show.legend = FALSE
    ) +
    geom_point(
      data = root_node,
      aes(
        x = .data$X,
        y = .data$Y
      ),
      shape = 21,
      size = 7.5,
      stroke = 1.1,
      color = "#2C7FB8",
      fill = "white"
    ) +
    scale_size_manual(
      values = count_sizes,
      breaks = count_levels,
      labels = count_labels,
      drop = FALSE,
      name = "Clone abundance"
    ) +
    scale_color_gradientn(
      colors = c(
        "#4682B4",
        "#4682B4",
        "#F4D8D4",
        "#C94A3A"
      ),
      values = scales::rescale(
        c(evar_limits[1], 0, 0.1, evar_limits[2]),
        from = evar_limits
      ),
      limits = evar_limits,
      oob = scales::squish,
      
      # E_change_var 缺失时也显示成蓝色，避免空白节点
      na.value = "#4682B4",
      
      name = "E_change_var",

    ) +
    coord_equal(
      xlim = c(
        -plotting_limit,
        plotting_limit
      ),
      ylim = c(
        -plotting_limit,
        plotting_limit
      ),
      clip = "off"
    ) +
    labs(
      title = paste0(
        group_name,
        "：14–21–28 d BCR clonal evolution"
      ),
      subtitle = paste0(
        "Radius = time; links = inferred minimum-Hamming lineage; ",
        "open ring = lineage not observed at a later time"
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 10,
        color = "#56666D",
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(
        t = 16,
        r = 16,
        b = 16,
        l = 16
      )
    ) +
    guides(
      size = guide_legend(
        order = 1,
        title.position = "top"
      ),
      color = guide_colorbar(
        order = 2,
        title.position = "top",
        barwidth = grid::unit(5, "cm")
      )
    )
  
  if (show_node_labels) {
    p <- p +
      geom_text(
        data = clone_nodes,
        aes(
          x = .data$X,
          y = .data$Y,
          label = .data$PatternIDs
        ),
        size = 2.2,
        nudge_y = 0.05,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }
  
  p
}

# -----------------------------
# 6. 单组与批量输出
# -----------------------------

build_one_group_temporal_tree <- function(nodes,
                                          group_name,
                                          days = c(14, 21, 28),
                                          max_hamming_distance = NULL,
                                          max_hamming_fraction = 0.10,
                                          evar_limits = c(-1, 1),
                                          show_node_labels = FALSE) {
  edges <- infer_temporal_edges(
    nodes = nodes,
    days = days,
    max_hamming_distance =
      max_hamming_distance,
    max_hamming_fraction =
      max_hamming_fraction
  )
  
  layout_nodes <- make_radial_time_layout(
    nodes = nodes,
    edges = edges,
    days = days
  )
  
  p <- plot_group_temporal_tree(
    group_name = group_name,
    layout_nodes = layout_nodes,
    edges = edges,
    days = days,
    evar_limits = evar_limits,
    show_node_labels = show_node_labels
  )
  
  list(
    nodes = nodes,
    layout_nodes = layout_nodes,
    edges = edges,
    plot = p,
    match_threshold = unique(
      edges$MatchThreshold
    )
  )
}

build_all_group_temporal_trees <- function(csv_dir,
                                           output_dir,
                                           days = c(14, 21, 28),
                                           top_n_per_time = 40,
                                           max_hamming_distance = NULL,
                                           max_hamming_fraction = 0.10,
                                           plot_width = 9,
                                           plot_height = 9,
                                           plot_dpi = 300,
                                           plot_device = c("png", "pdf"),
                                           show_node_labels = FALSE,
                                           encoding = "UTF-8",
                                           strict_filename_parsing = FALSE) {
  plot_device <- match.arg(plot_device)
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  manifest <- build_file_manifest(
    csv_dir = csv_dir,
    days = days,
    strict_filename_parsing =
      strict_filename_parsing
  )
  
  readr::write_csv(
    manifest,
    file.path(
      output_dir,
      "BCR_14_21_28_file_manifest.csv"
    )
  )
  
  message("成功解析的文件：")
  print(
    manifest %>%
      select(
        .data$Sample,
        .data$Group,
        .data$Day
      )
  )
  
  all_data <- read_manifest_data(
    manifest,
    encoding = encoding
  )
  
  group_names <- sort(
    unique(all_data$Group)
  )
  
  valid_groups <- character(0)
  skipped_groups <- character(0)
  
  for (group_name in group_names) {
    group_days <- sort(
      unique(
        all_data$Day[
          all_data$Group == group_name
        ]
      )
    )
    
    if (all(days %in% group_days)) {
      valid_groups <- c(
        valid_groups,
        group_name
      )
    } else {
      skipped_groups <- c(
        skipped_groups,
        group_name
      )
      warning(
        "跳过组 ", group_name,
        "：缺少时间点 ",
        paste(
          setdiff(days, group_days),
          collapse = ", "
        )
      )
    }
  }
  
  if (length(valid_groups) == 0) {
    stop(
      "没有任何组同时包含 ",
      paste(days, collapse = "/"),
      " 天数据。"
    )
  }
  
  prepared_nodes <- lapply(
    valid_groups,
    function(group_name) {
      prepare_group_nodes(
        group_df = all_data %>%
          filter(
            .data$Group == group_name
          ),
        group_name = group_name,
        days = days,
        top_n_per_time = top_n_per_time
      )
    }
  )
  names(prepared_nodes) <- valid_groups
  
  global_evar_limits <- get_global_evar_limits(
    prepared_nodes
  )
  
  message(
    "所有组统一 E_change_var 色阶：",
    paste(
      signif(global_evar_limits, 4),
      collapse = " ~ "
    )
  )
  
  results <- vector(
    "list",
    length(valid_groups)
  )
  names(results) <- valid_groups
  
  for (group_name in valid_groups) {
    message(
      "正在构建三时间点联合谱系树：",
      group_name
    )
    
    result <- build_one_group_temporal_tree(
      nodes = prepared_nodes[[group_name]],
      group_name = group_name,
      days = days,
      max_hamming_distance =
        max_hamming_distance,
      max_hamming_fraction =
        max_hamming_fraction,
      evar_limits = global_evar_limits,
      show_node_labels = show_node_labels
    )
    
    safe_group_name <- gsub(
      "[^\\p{L}\\p{N}_-]+",
      "_",
      group_name,
      perl = TRUE
    )
    
    plot_file <- file.path(
      output_dir,
      paste0(
        safe_group_name,
        "_BCR_14_21_28_lineage_tree.",
        plot_device
      )
    )
    
    ggsave(
      filename = plot_file,
      plot = result$plot,
      width = plot_width,
      height = plot_height,
      dpi = plot_dpi,
      bg = "white"
    )
    
    readr::write_csv(
      result$layout_nodes,
      file.path(
        output_dir,
        paste0(
          safe_group_name,
          "_BCR_14_21_28_nodes.csv"
        )
      )
    )
    
    readr::write_csv(
      result$edges,
      file.path(
        output_dir,
        paste0(
          safe_group_name,
          "_BCR_14_21_28_edges.csv"
        )
      )
    )
    
    message("树图已保存：", plot_file)
    message(
      "本组允许的最大 Hamming 距离：",
      result$match_threshold
    )
    
    results[[group_name]] <- result
  }
  
  attr(results, "manifest") <- manifest
  attr(results, "global_evar_limits") <-
    global_evar_limits
  attr(results, "skipped_groups") <-
    skipped_groups
  
  results
}

# -----------------------------
# 7. 实际运行区
# -----------------------------

csv_dir <- "Your GC folder"
out_dir <- file.path(
  csv_dir,
  "GC_BCR_lineage_trees_30"
)

bcr_temporal_tree_results <-
  build_all_group_temporal_trees(
    csv_dir = csv_dir,
    output_dir = out_dir,
    days = c(14, 21, 28),
    
    # 每个时间点最多保留多少个 clone。
    top_n_per_time = 30,
    
    # 相似谱系判定：
    # 默认允许 char 总长度的 10% 发生变化。
    # 例如 char 长度为 50，则最大 Hamming 距离为 5。
    max_hamming_fraction = 0.1,
    
    # 如果希望直接指定最多允许几个不同位点，
    # 可把下一行改成 2、3 等；NULL 表示使用上面的比例。
    max_hamming_distance = NULL,
    
    plot_width = 9,
    plot_height = 9,
    plot_dpi = 300,
    plot_device = "png",
    show_node_labels = FALSE,
    strict_filename_parsing = FALSE
  )

