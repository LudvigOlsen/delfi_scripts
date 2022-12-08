
library(tidyverse)
library(caret)
library(pROC)
library(optparse)
library(stringr)

# Command line arguments
option_list = list(
  make_option(
    c("-i", "--in_features_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A file with the short+total features. Must have the extension `.csv`."
  ),
  make_option(
    c("-m", "--in_meta_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A file where the first 2 columns contain the sample IDs and labels. Must have the extension `.csv`."
  ),
  make_option(
    c("-o", "--out_model_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The filepath to save the model at. Must have the extension `.rds`."
  ),
  make_option(
    c("-p", "--out_preds_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "The filepath to save the predictions at. Must have the extension `.csv`."
  ),
  make_option(
    c("-c", "--control_labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = "String with comma-separated control labels. Replace whitespace with underscores."
  ),
  make_option(
    c("-l", "--cancer_labels"),
    action = "store",
    default = NA,
    type = 'character',
    help = "String with comma-separated cancer labels. Replace whitespace with underscores."
  )
)
opt = parse_args(OptionParser(option_list = option_list))
if (is.na(opt$in_features_file)){
  stop("--in_features_file was not specified.")
}
if (is.na(opt$in_meta_file)){
  stop("--in_meta_file was not specified.")
}
if (is.na(opt$out_model_file)){
  stop("--out_model_file was not specified.")
}
if (is.na(opt$out_preds_file)){
  stop("--out_preds_file was not specified.")
}
if (is.na(opt$control_labels)){
  stop("--control_labels was not specified.")
}
if (is.na(opt$cancer_labels)){
  stop("--cancer_labels was not specified.")
}
if (str_sub(opt$in_features_file, start= -4) != ".csv"){
  stop("--in_features_file must have the extension '.csv'.")
}
if (str_sub(opt$in_meta_file, start= -4) != ".csv"){
  stop("--in_meta_file must have the extension '.csv'.")
}
if (str_sub(opt$out_model_file, start= -4) != ".rds"){
  stop("--out_model_file must have the extension '.rds'.")
}
if (str_sub(opt$out_preds_file, start= -4) != ".csv"){
  stop("--out_preds_file must have the extension '.csv'.")
}

# Load in features and add labels
features.sl <- read_csv(opt$in_features_file)
meta_data <- read_csv(opt$in_meta_file)
# Assign class to data frame (but without spaces)
original_labels <- stringr::str_replace_all(meta_data[[2]], " ", "_")
control_labels <- str_split(opt$control_labels, pattern = ",")
cancer_labels <- str_split(opt$cancer_labels, pattern = ",")

print(paste0("Control labels: ", control_labels, collapse = TRUE))
print(paste0("Cancer labels: ", cancer_labels, collapse = TRUE))

features.sl["type"] <- dplyr::case_when(
  original_labels %in% control_labels ~ "Control",
  original_labels %in% cancer_labels ~ "Cancer",
  TRUE ~ "Discard"
)

features.sl <- features.sl %>%
  dplyr::filter(type != "Discard")

print("Loaded features and labels: ")
print(head(features.sl, 10))

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  verboseIter = FALSE,
  savePredictions = TRUE,
  classProbs = TRUE,
  # preProcOptions=list(thres = 0.90),
  summaryFunction = twoClassSummary
)

# ####### Only short/total coverage
set.seed(1234)
model_sl <- caret::train(
  type ~ .,
  data = features.sl,
  method = 'gbm',
  tuneGrid = data.frame(
    n.trees = 150,
    interaction.depth = 3,
    shrinkage = 0.1,
    n.minobsinnode = 10
  ),
  preProcess = c("corr", "nzv"),
  trControl = ctrl
)

# Save models
models.list <- list("SL"=model_sl) #"all"=model_gbm, "z"=model_z)
saveRDS(models.list, opt$out_model_file)

pred.tbl <- model_sl$pred %>%
  dplyr::as_tibble() %>%
  filter(n.trees==150, interaction.depth==3)

# Save predictions before averaging them
write.csv(pred.tbl, file = opt$out_preds_file)

# Average predictions and add sample
pred.tbl <- pred.tbl %>%
  group_by(rowIndex) %>%
  dplyr::summarize(
    obs=obs[1],
    Cancer=mean(Cancer)
  )
pred.tbl$sample <-  meta_data[[1]]

# Save averaged predictions
write.csv(pred.tbl, file = paste0(str_sub(opt$out_preds_file, end= -5)), ".avg.csv")

# ## 95% specificity
# cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
#            arrange(desc(Cancer)))$Cancer[11]
# cutoff98 <- (pred.tbl %>% filter(type=="Healthy") %>%
#              arrange(desc(Cancer)))$Cancer[5]
# ## 90% specificity cutoff to be used in tissue prediction.
# cutoff90 <- (pred.tbl %>% filter(type=="Healthy") %>%
#              arrange(desc(Cancer)))$Cancer[21]
#
# pred.tbl <- pred.tbl %>%
#     mutate(detected95 = ifelse(Cancer > cutoff, "Detected", "Not detected"),
#        detected98 = ifelse(Cancer > cutoff98, "Detected", "Not detected"),
#        detected90 = ifelse(Cancer > cutoff90, "Detected", "Not detected"),
#        stage = gsub("A|B|C", "", `Stage at Diagnosis`))
#
# write.csv(inner_join(summary.df %>% select(-contains("Z Score")), pred.tbl %>%
#                      select(rowIndex, sample, stage, Cancer, detected95, detected98),
#                      by=c("sample"="sample")),"../inst/extdata/predictions_gbm.csv",
#           row.names=FALSE)
