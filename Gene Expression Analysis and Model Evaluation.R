############################################################
# Load Packages
############################################################
library(rtracklayer)
library(GenomicRanges)
library(stats)
library(BiocParallel)
library(glmnet)
library(randomForest)
library(earth)
library(pROC)
library(e1071)
library(ggplot2)

############################################################
# Parallel Setup
############################################################
num_workers <- 18
register(SnowParam(workers = num_workers))

############################################################
# Files
############################################################
histone_files <- c(
  "wgEncodeBroadHistoneGm12878H3k4me3StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k9acStdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k36me3StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H2azStdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k79me2StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k4me1StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k4me2StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k9me3StdAln_3Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H4k20me1StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k27me3StdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeBroadHistoneGm12878H3k27acStdAln_2Reps.norm5.rawsignal.bw",
  "wgEncodeUwDnaseGm12878Aln_2Reps.norm5.rawsignal.bw"
)

gtf_file <- "gencode.v7.annotation_updated_ncrna_host.gtf"
cage_bw_file <- "wgEncodeRikenCageGm12878NucleusPapPlusSignalRep1.bigWig"

############################################################
# Parameters
############################################################
upstream <- 2000
downstream <- 2000
nbins_up <- 40
nbins_down <- 40
window_size <- 50
min_length <- 4100
fixed_pseudocount <- 0.1

############################################################
# Data Processing (Gene-Level)
############################################################
gtf_data <- import(gtf_file, format = "gtf")
gtf_data <- gtf_data[seqnames(gtf_data) == "chr1"]

transcripts <- gtf_data[gtf_data$type == "transcript"]
transcript_lengths <- width(transcripts)
valid_transcripts <- transcripts[transcript_lengths >= min_length]

tss_positions <- ifelse(strand(valid_transcripts) == "+", start(valid_transcripts), end(valid_transcripts))
tss_gr <- GRanges(seqnames(valid_transcripts), IRanges(tss_positions, width=1), strand=strand(valid_transcripts), gene_id=valid_transcripts$gene_id)
tss_gr <- unique(tss_gr)

cage_data <- import(cage_bw_file, format = "bigWig")
cage_data <- cage_data[seqnames(cage_data) == "chr1"]

promoters <- function(x, upstream=2000, downstream=200) {
  out <- x
  plus_idx <- which(strand(x) == "+")
  minus_idx <- which(strand(x) == "-")
  start(out)[plus_idx] <- start(x)[plus_idx] - upstream
  end(out)[plus_idx]   <- start(x)[plus_idx] + downstream - 1
  start(out)[minus_idx] <- end(x)[minus_idx] - downstream + 1
  end(out)[minus_idx]   <- end(x)[minus_idx] + upstream
  trim(out)
}

tss_windows <- promoters(tss_gr, upstream=window_size, downstream=window_size+1)
ov_cage <- findOverlaps(tss_windows, cage_data)
cage_sum <- tapply(cage_data$score[subjectHits(ov_cage)], queryHits(ov_cage), sum, default=0)
cage_signal <- numeric(length(tss_windows))
cage_signal[as.integer(names(cage_sum))] <- cage_sum
total_cage_tags <- sum(cage_data$score)
cage_signal <- cage_signal / (total_cage_tags / 1e6)

gene_tss_list <- split(tss_gr, tss_gr$gene_id)
gene_signals_list <- split(cage_signal, tss_gr$gene_id)
max_expr_idx_per_gene <- sapply(gene_signals_list, which.max)
selected_tss_list <- mapply(function(tss_grp, idx) tss_grp[idx], tss_grp = gene_tss_list, idx = max_expr_idx_per_gene, SIMPLIFY = FALSE)
selected_signal <- mapply(function(sig_grp, idx) sig_grp[idx], sig_grp = gene_signals_list, idx = max_expr_idx_per_gene)
selected_tss_list_gr <- GRangesList(selected_tss_list)
selected_tss <- unlist(selected_tss_list_gr)

create_bins <- function(tss, upstream=2000, downstream=2000, nbins_up=40, nbins_down=40) {
  strand_tss <- as.character(strand(tss))
  pos <- start(tss)
  if (strand_tss == "-") {
    region_start <- pos - downstream
    region_end   <- pos + upstream
  } else {
    region_start <- pos - upstream
    region_end   <- pos + downstream
  }
  region <- GRanges(seqnames(tss), IRanges(min(region_start, region_end), max(region_start, region_end)), strand=strand_tss)
  total_bins <- nbins_up + nbins_down + 1
  bin_breaks <- seq(start(region), end(region), length.out = total_bins + 1)
  bins <- GRanges(seqnames(region), IRanges(start=bin_breaks[-(total_bins+1)], end=bin_breaks[-1]), strand=strand_tss)
  trim(bins)
}

all_bins_list <- lapply(seq_along(selected_tss), function(i) {
  create_bins(selected_tss[i], upstream, downstream, nbins_up, nbins_down)
})

set.seed(123)
n_genes <- length(selected_tss)
D1_size <- round(n_genes / 3)
D1_idx <- sample(seq_len(n_genes), size = D1_size)
D2_idx <- setdiff(seq_len(n_genes), D1_idx)

D1_bins_list <- all_bins_list[D1_idx]
D1_cage_signal <- selected_signal[D1_idx]

D2_bins_list <- all_bins_list[D2_idx]
D2_cage_signal <- selected_signal[D2_idx]

mean_bin_signal <- function(bins, cov) {
  bins_by_chr <- split(bins, seqnames(bins))
  result <- numeric(length(bins))
  start_idx <- 1
  for (chr in names(bins_by_chr)) {
    chr_bins <- bins_by_chr[[chr]]
    if (!chr %in% names(cov)) {
      end_idx <- start_idx + length(chr_bins) - 1
      result[start_idx:end_idx] <- NA
      start_idx <- end_idx + 1
      next
    }
    chr_cov <- cov[[chr]]
    chr_bins <- trim(chr_bins)
    if (length(chr_bins) == 0) next
    v <- Views(chr_cov, start=start(chr_bins), end=end(chr_bins))
    means <- viewMeans(v)
    if (length(means) == 0) means <- rep(NA, length(chr_bins))
    end_idx <- start_idx + length(chr_bins) - 1
    result[start_idx:end_idx] <- means
    start_idx <- end_idx + 1
  }
  result
}

find_best_bin_fixed_pc <- function(hist_cov, bins_list, expression, pc=0.1) {
  signals_list <- lapply(bins_list, mean_bin_signal, cov=hist_cov)
  if (length(signals_list) == 0) return(list(pseudocount=pc, best_bin=1))
  signal_matrix <- do.call(rbind, signals_list)
  if (is.null(signal_matrix) || nrow(signal_matrix) == 0) return(list(pseudocount=pc, best_bin=1))
  if (all(is.na(signal_matrix))) return(list(pseudocount=pc, best_bin=1))
  log_mat <- log2(signal_matrix + pc)
  cors <- apply(log_mat, 2, function(x) cor(x, expression, use="pairwise.complete.obs"))
  if (all(is.na(cors))) return(list(pseudocount=pc, best_bin=1))
  current_best_bin_idx <- which.max(abs(cors))
  if (length(current_best_bin_idx)==0 || is.na(current_best_bin_idx)) current_best_bin_idx <- 1
  list(pseudocount=pc, best_bin=current_best_bin_idx)
}

hist_cov_list_d1 <- lapply(histone_files, function(hfile) {
  hist_data <- import(hfile, format="bigWig")
  hist_data <- hist_data[seqnames(hist_data) == "chr1"]
  coverage(hist_data, weight="score")
})

pc_and_bin_list <- list()
for (i in seq_along(histone_files)) {
  hfile <- histone_files[i]
  hist_cov <- hist_cov_list_d1[[i]]
  res <- find_best_bin_fixed_pc(hist_cov, D1_bins_list, D1_cage_signal, pc=fixed_pseudocount)
  pc_and_bin_list[[hfile]] <- res
}

############################################################
# Prepare final dataset (D2) for CV
############################################################
hist_cov_list_d2 <- lapply(histone_files, function(hfile) {
  hist_data <- import(hfile, format="bigWig")
  hist_data <- hist_data[seqnames(hist_data) == "chr1"]
  coverage(hist_data, weight="score")
})

make_feature_matrix <- function(histone_files, bins_list, cage_signal, pc_and_bin, hist_cov_list) {
  signals_list <- lapply(seq_along(histone_files), function(i) {
    hfile <- histone_files[i]
    cov <- hist_cov_list[[i]]
    mat_list <- lapply(bins_list, mean_bin_signal, cov=cov)
    mat <- do.call(rbind, mat_list)
    pc <- pc_and_bin[[hfile]]$pseudocount
    bbin <- pc_and_bin[[hfile]]$best_bin
    log_mat <- log2(mat + pc)
    log_mat[, bbin, drop=FALSE]
  })
  X <- do.call(cbind, signals_list)
  colnames(X) <- histone_files
  X
}

X_all <- make_feature_matrix(histone_files, D2_bins_list, D2_cage_signal, pc_and_bin_list, hist_cov_list_d2)
Y_all <- rep(NA, length(D2_cage_signal))
nonzero_idx <- which(D2_cage_signal > 0)
Y_all[nonzero_idx] <- log2(D2_cage_signal[nonzero_idx])
Y_all[is.na(Y_all)] <- 0
is_on_all <- Y_all > 0.1

############################################################
# 10-fold CV
############################################################
set.seed(123)
k <- 10
n <- nrow(X_all)
folds <- sample(rep(1:k, length.out = n))

misclass_rate <- function(true_labels, pred_labels) mean(true_labels != pred_labels)
rmse_log <- function(true_log, pred_log) sqrt(mean((true_log - pred_log)^2, na.rm=TRUE))
pearson_r <- function(true_log, pred_log) cor(true_log, pred_log, use="complete.obs")

two_step_combos <- c("Logistic_Lasso","Logistic_RF_Reg","Logistic_M","Logistic_Sv",
                     "RF_Clf_Lasso","RF_Clf_RF_Reg","RF_Clf_M","RF_Clf_Sv",
                     "SVM_Clf_Lasso","SVM_Clf_RF_Reg","SVM_Clf_M","SVM_Clf_Sv")

auc_values <- matrix(0, nrow=k, ncol=3, dimnames=list(NULL,c("Logistic","RF_Clf","SVM_Clf")))
misclass_values <- matrix(0, nrow=k, ncol=3, dimnames=list(NULL,c("Logistic","RF_Clf","SVM_Clf")))
rmse_values <- matrix(0, nrow=k, ncol=length(two_step_combos), dimnames=list(NULL,two_step_combos))
r_values <- matrix(0, nrow=k, ncol=length(two_step_combos), dimnames=list(NULL,two_step_combos))

rf_roc_for_plot <- NULL
chosen_combo <- "RF_Clf_RF_Reg"
chosen_iteration <- 1
store_actual <- NULL
store_pred <- NULL

for(i in 1:k) {
  test_idx <- which(folds == i)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  X_train_cv <- X_all[train_idx, ]
  Y_train_cv <- Y_all[train_idx]
  is_on_train_cv <- is_on_all[train_idx]
  
  X_test_cv <- X_all[test_idx, ]
  Y_test_cv <- Y_all[test_idx]
  is_on_test_cv <- is_on_all[test_idx]
  
  clf_train_data_cv <- data.frame(X_train_cv, is_on = factor(is_on_train_cv, levels=c(FALSE,TRUE)))
  clf_test_data_cv  <- data.frame(X_test_cv, is_on = factor(is_on_test_cv, levels=c(FALSE,TRUE)))
  
  # Classification
  log_model_cv <- glm(is_on ~ ., data=clf_train_data_cv, family="binomial")
  log_pred_prob_cv <- predict(log_model_cv, newdata=clf_test_data_cv, type="response")
  log_pred_class_cv <- log_pred_prob_cv > 0.5
  
  rf_clf_cv <- randomForest(is_on ~ ., data=clf_train_data_cv)
  rf_clf_prob_cv <- predict(rf_clf_cv, newdata=clf_test_data_cv, type="prob")[,2]
  rf_clf_pred_cv <- rf_clf_prob_cv > 0.5
  rf_roc_cv <- roc(is_on_test_cv, rf_clf_prob_cv)
  if (i == 1) { rf_roc_for_plot <- rf_roc_cv }
  
  svm_clf_cv <- svm(is_on ~ ., data=clf_train_data_cv, probability=TRUE, kernel="radial")
  svm_clf_pred_prob_cv <- attr(predict(svm_clf_cv, newdata=clf_test_data_cv, probability=TRUE),"probabilities")[,2]
  svm_clf_pred_cv <- svm_clf_pred_prob_cv > 0.5
  
  log_roc_cv <- roc(is_on_test_cv, log_pred_prob_cv)
  auc_values[i,"Logistic"] <- auc(log_roc_cv)
  misclass_values[i,"Logistic"] <- misclass_rate(is_on_test_cv, log_pred_class_cv)
  
  rf_roc_cv <- roc(is_on_test_cv, rf_clf_prob_cv)
  auc_values[i,"RF_Clf"] <- auc(rf_roc_cv)
  misclass_values[i,"RF_Clf"] <- misclass_rate(is_on_test_cv, rf_clf_pred_cv)
  
  svm_roc_cv <- roc(is_on_test_cv, svm_clf_pred_prob_cv)
  auc_values[i,"SVM_Clf"] <- auc(svm_roc_cv)
  misclass_values[i,"SVM_Clf"] <- misclass_rate(is_on_test_cv, svm_clf_pred_cv)
  
  on_train_cv_idx <- which(is_on_train_cv)
  on_test_cv_idx <- which(is_on_test_cv)
  
  Y_train_on_cv <- Y_train_cv[on_train_cv_idx]
  Y_test_on_cv <- Y_test_cv[on_test_cv_idx]
  X_train_on_cv <- X_train_cv[on_train_cv_idx,]
  X_test_on_cv <- X_test_cv[on_test_cv_idx,]
  
  lasso_model_cv <- cv.glmnet(X_train_on_cv, Y_train_on_cv, alpha=1, family="gaussian")
  best_lambda_cv <- lasso_model_cv$lambda.min
  Y_pred_lasso_cv <- predict(lasso_model_cv, newx=X_test_on_cv, s=best_lambda_cv)
  
  rf_reg_cv <- randomForest(x=X_train_on_cv, y=Y_train_on_cv)
  Y_pred_rf_reg_cv <- predict(rf_reg_cv, newdata=X_test_on_cv)
  
  mars_reg_cv <- earth(Y_train_on_cv ~ ., data=data.frame(X_train_on_cv, Y_train_on_cv=Y_train_on_cv), degree=1)
  Y_pred_mars_cv <- predict(mars_reg_cv, newdata=data.frame(X_test_on_cv))
  
  svm_reg_cv <- svm(Y_train_on_cv ~ ., data=data.frame(X_train_on_cv, Y_train_on_cv=Y_train_on_cv), kernel="radial")
  Y_pred_svm_reg_cv <- predict(svm_reg_cv, newdata=data.frame(X_test_on_cv))
  
  combos <- list(
    Logistic_Lasso  = list(on_idx = which(log_pred_class_cv), pred=Y_pred_lasso_cv),
    Logistic_RF_Reg = list(on_idx = which(log_pred_class_cv), pred=Y_pred_rf_reg_cv),
    Logistic_M      = list(on_idx = which(log_pred_class_cv), pred=Y_pred_mars_cv),
    Logistic_Sv     = list(on_idx = which(log_pred_class_cv), pred=Y_pred_svm_reg_cv),
    RF_Clf_Lasso    = list(on_idx = which(rf_clf_pred_cv), pred=Y_pred_lasso_cv),
    RF_Clf_RF_Reg   = list(on_idx = which(rf_clf_pred_cv), pred=Y_pred_rf_reg_cv),
    RF_Clf_M        = list(on_idx = which(rf_clf_pred_cv), pred=Y_pred_mars_cv),
    RF_Clf_Sv       = list(on_idx = which(rf_clf_pred_cv), pred=Y_pred_svm_reg_cv),
    SVM_Clf_Lasso   = list(on_idx = which(svm_clf_pred_cv), pred=Y_pred_lasso_cv),
    SVM_Clf_RF_Reg  = list(on_idx = which(svm_clf_pred_cv), pred=Y_pred_rf_reg_cv),
    SVM_Clf_M       = list(on_idx = which(svm_clf_pred_cv), pred=Y_pred_mars_cv),
    SVM_Clf_Sv      = list(on_idx = which(svm_clf_pred_cv), pred=Y_pred_svm_reg_cv)
  )
  
  for(combo_name in names(combos)) {
    mod_cv <- combos[[combo_name]]
    on_predicted_on_idx <- intersect(on_test_cv_idx, mod_cv$on_idx)
    rel_idx <- match(on_predicted_on_idx, on_test_cv_idx)
    rel_idx <- rel_idx[!is.na(rel_idx)]
    if(length(rel_idx) > 0) {
      rmse_values[i, combo_name] <- rmse_log(Y_test_on_cv[rel_idx], mod_cv$pred[rel_idx])
      r_values[i, combo_name] <- pearson_r(Y_test_on_cv[rel_idx], mod_cv$pred[rel_idx])
      if (combo_name == chosen_combo && i == chosen_iteration) {
        store_actual <- Y_test_on_cv[rel_idx]
        store_pred <- mod_cv$pred[rel_idx]
      }
    } else {
      rmse_values[i, combo_name] <- NA
      r_values[i, combo_name] <- NA
    }
  }
}

mean_auc <- colMeans(auc_values)
mean_misclass <- colMeans(misclass_values, na.rm=TRUE)
mean_rmse <- colMeans(rmse_values, na.rm=TRUE)
mean_r <- colMeans(r_values, na.rm=TRUE)

cat("Classification Comparison (Mean over 10 folds):\n")
cat("Model       Misclass   AUC\n")
cat("Logistic    ", mean_misclass["Logistic"], "  ", mean_auc["Logistic"], "\n")
cat("RF_Clf      ", mean_misclass["RF_Clf"], "  ", mean_auc["RF_Clf"], "\n")
cat("SVM_Clf     ", mean_misclass["SVM_Clf"], "  ", mean_auc["SVM_Clf"], "\n")

cat("\nRegression Comparison (Mean over 10 folds):\n")
cat("Model                      RMSE_log   Pearson_r\n")
for(combo_name in two_step_combos) {
  cat(combo_name, "  ", mean_rmse[combo_name], "   ", mean_r[combo_name], "\n")
}

if(!is.null(rf_roc_for_plot)) {
  plot(rf_roc_for_plot,
       main="Random Forest ROC Curve (Fold 1)",
       legacy.axes=TRUE,
       xlab="False Positive Rate",
       ylab="True Positive Rate")
}

if(!is.null(store_actual) && !is.null(store_pred)) {
  corr <- cor(store_actual, store_pred, use="complete.obs")
  df_plot <- data.frame(Actual=store_actual, Predicted=store_pred)
  ggplot(df_plot, aes(x=Actual, y=Predicted)) +
    geom_point(color="blue") +
    geom_smooth(method="lm", color="red", linetype="dashed", se=FALSE) +
    labs(title=paste("Predicted vs. Actual for", chosen_combo, "(Fold", chosen_iteration,")"),
         subtitle=paste("Pearson's r =", round(corr, 3)),
         x="Actual",
         y="Predicted"
    )
}

############################################################
# Variable Importance
############################################################
seeds <- c(111, 222, 333, 444, 555, 666, 777, 888, 999, 1010)
clf_importances <- list()
reg_importances <- list()
on_idx_full <- which(is_on_all)
X_on_full <- X_all[on_idx_full,]
Y_on_full <- Y_all[on_idx_full]

for(i in seq_along(seeds)) {
  set.seed(seeds[i])
  clf_full_data <- data.frame(X_all, is_on = factor(is_on_all, levels=c(FALSE,TRUE)))
  rf_clf_full <- randomForest(is_on ~ ., data=clf_full_data)
  clf_imp <- importance(rf_clf_full)
  clf_mdg <- clf_imp[,"MeanDecreaseGini"]
  clf_importances[[i]] <- clf_mdg
  
  rf_reg_full <- randomForest(x=X_on_full, y=Y_on_full)
  reg_imp <- importance(rf_reg_full)
  reg_inp <- reg_imp[,"IncNodePurity"]
  reg_importances[[i]] <- reg_inp
}

clf_mat <- do.call(cbind, clf_importances)
reg_mat <- do.call(cbind, reg_importances)

mean_clf_importance <- rowMeans(clf_mat, na.rm=TRUE)
mean_reg_importance <- rowMeans(reg_mat, na.rm=TRUE)

cat("\nAveraged RF Classification MeanDecreaseGini:\n")
clf_var_imp_df <- data.frame(Feature=rownames(clf_imp), MeanDecreaseGini=mean_clf_importance)
clf_var_imp_df <- clf_var_imp_df[order(clf_var_imp_df$MeanDecreaseGini, decreasing=TRUE), ]
print(clf_var_imp_df)

cat("\nAveraged RF Regression IncNodePurity:\n")
reg_var_imp_df <- data.frame(Feature=rownames(reg_imp), IncNodePurity=mean_reg_importance)
reg_var_imp_df <- reg_var_imp_df[order(reg_var_imp_df$IncNodePurity, decreasing=TRUE), ]
print(reg_var_imp_df)

extract_short_name <- function(fname) {
  if (grepl("UwDnase", fname)) {
    return("Dnase")
  } else {
    part <- sub(".*Gm12878", "", fname)
    part <- sub("StdAln.*", "", part)
    return(part)
  }
}

clf_features <- rownames(clf_mat)
num_features_clf <- length(clf_features)
num_runs_clf <- ncol(clf_mat)
clf_df_list <- list()
for (r in 1:num_runs_clf) {
  run_df <- data.frame(Feature = clf_features, Importance = clf_mat[,r], Run = factor(r))
  clf_df_list[[r]] <- run_df
}
clf_long_df <- do.call(rbind, clf_df_list)
clf_long_df$Feature <- sapply(clf_long_df$Feature, extract_short_name)
mean_clf_imp <- tapply(clf_long_df$Importance, clf_long_df$Feature, mean, na.rm=TRUE)
order_feats_clf <- names(sort(mean_clf_imp, decreasing=TRUE))
clf_long_df$Feature <- factor(clf_long_df$Feature, levels=order_feats_clf)

reg_features <- rownames(reg_mat)
num_features_reg <- length(reg_features)
num_runs_reg <- ncol(reg_mat)
reg_df_list <- list()
for (r in 1:num_runs_reg) {
  run_df <- data.frame(Feature = reg_features, Importance = reg_mat[,r], Run = factor(r))
  reg_df_list[[r]] <- run_df
}
reg_long_df <- do.call(rbind, reg_df_list)
reg_long_df$Feature <- sapply(reg_long_df$Feature, extract_short_name)
mean_reg_imp <- tapply(reg_long_df$Importance, reg_long_df$Feature, mean, na.rm=TRUE)
order_feats_reg <- names(sort(mean_reg_imp, decreasing=TRUE))
reg_long_df$Feature <- factor(reg_long_df$Feature, levels=order_feats_reg)

ggplot(clf_long_df, aes(x=Feature, y=Importance)) +
  geom_boxplot(fill="lightblue") +
  labs(title="Variable Importance for classification",
       x="Feature",
       y="MeanDecreaseGini") +
  theme(axis.text.x = element_text(angle=90, hjust=1))

ggplot(reg_long_df, aes(x=Feature, y=Importance)) +
  geom_boxplot(fill="lightgreen") +
  labs(title="Variable Importance for Regression",
       x="Feature",
       y="IncNodePurity") +
  theme(axis.text.x = element_text(angle=90, hjust=1))