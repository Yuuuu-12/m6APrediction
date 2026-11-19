# Functions reused from Practical 5. Add documentation in your new package.
#' Encode DNA Sequence Features
#'
#' Encodes the input DNA sequences into a numeric representation suitable for the random forest model.
#'
#' @param dna_strings The input DNA sequence(s).
#' @return A numeric matrix or vector representing the encoded sequence.
#' @importFrom stats predict
#' @examples
#' # No examples for internal helper function

dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Sequences
#'
#' Predicts m6A probability and status for a batch of sequences.
#'
#' @param ml_fit The fitted random forest model object.
#' @param feature_df Dataframe or matrix containing the input features.
#' @param positive_threshold Numeric threshold for classification. Scores strictly greater than the threshold are Positive; equality is Negative.
#' @return A dataframe with added columns `predicted_m6A_prob` and `predicted_m6A_status`.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' # Load the fitted model and example data
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction", mustWork = TRUE))
#'
#' # Run batch prediction on the first 3 rows to keep output short
#' batch_check <- prediction_multiple(ml_fit, example_df[1:3, ], positive_threshold = 0.6)
#' print(batch_check[, c("predicted_m6A_prob", "predicted_m6A_status")])



prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A Site for a Single Sequence
#'
#' Predicts m6A status and returns a named character vector.
#'
#' @param ml_fit The fitted random forest model object.
#' @param gc_content Numeric, GC content.
#' @param RNA_type Character, RNA type (e.g., "mRNA").
#' @param RNA_region Character, RNA region (e.g., "CDS").
#' @param exon_length Numeric, exon length.
#' @param distance_to_junction Numeric, distance to junction.
#' @param evolutionary_conservation Numeric, evolutionary conservation score.
#' @param DNA_5mer Character, 5-mer DNA sequence (e.g., "GGACA").
#' @param positive_threshold Numeric threshold for classification.
#' @return A named character vector containing `predicted_m6A_prob` and `predicted_m6A_status`.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' # Load the fitted model
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE))
#'
#' # Predict a single 5-base context sequence (GGACA)
#' prediction_single(ml_fit,
#'                   gc_content = 0.6,
#'                   RNA_type = "mRNA",
#'                   RNA_region = "CDS",
#'                   exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)


prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
