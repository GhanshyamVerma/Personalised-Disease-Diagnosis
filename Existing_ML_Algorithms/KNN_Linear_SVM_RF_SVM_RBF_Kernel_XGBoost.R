# Required Programming Language
# R version 3.6.2 or above

# Required packages
required_packages <- c("caret", "dplyr", "e1071", "ggplot2", "class", 
                       "randomForest", "pROC", "kernlab", "xgboost")

# Check if required r libraries are already installed or not
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# If required r packages or libraries are not already installed, Install them.
if(length(new_packages)) install.packages(new_packages)

# load packages# Required Programming Language
# R version 3.6.2 or above

# Required packages
required_packages <- c("caret", "dplyr", "e1071", "ggplot2", "class", 
                       "randomForest", "pROC", "kernlab", "xgboost")

# Check if required r libraries are already installed or not
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# If required r packages or libraries are not already installed, Install them.
if(length(new_packages)) install.packages(new_packages)

# load packages
library(caret) 
library(dplyr) 
library(e1071) 
library(ggplot2)
library(class) 
library(randomForest) 
library(pROC) 
library(kernlab) 
library(xgboost)

# Provide seed values for data partition and ml model building

data_partition_seed <- 7
ml_model_seed <- 1234


# Define an input path
input_path <- "./Datasets/" 

# Define an output path
output_path <- "./All_ML_Models_Results/" 

# Read Gene Expression Datasets
read_gene_expression_data <- function(filename,input_path) {
  path_filename <- paste0(input_path, filename)
  print(paste0("Reading CSV file called ", path_filename))
  return(read.csv(file = path_filename, header = TRUE, sep=",", check.names = FALSE))
}


# Split data into train, validation, and test
split_data <- function(Gene_Exp_Data, data_partition_seed, dataset_name) {
  total_data <- Gene_Exp_Data
  Data_target <- total_data %>% filter(Time > 0)
  
  # Split into train and test
  set.seed(data_partition_seed)
  index_test <- createDataPartition(y = Data_target$True_Class_Label,
                                    p = 0.50, list = FALSE)
  test_data <- Data_target[index_test, ]
  train_data <- Data_target[-index_test, ]
  
  dim(train_data)
  dim(test_data)
  
  # Split test data into holdout test and validation
  set.seed(data_partition_seed)
  index_test <- createDataPartition(y = test_data$True_Class_Label,
                                    p = 0.50, list = FALSE)
  hold_out_test <- test_data[index_test, ]
  valid_data <- test_data[-index_test, ]
  dim(hold_out_test)
  dim(valid_data)
  
  # Full train data for final model building
  full_train_data <- rbind(train_data,valid_data)
  
  # train, valid and holdout_test all time points
  g_train_data <- total_data %>% filter(Super_Subject_ID %in% train_data$Super_Subject_ID)
  g_valid_data <- total_data %>% filter(Super_Subject_ID %in% valid_data$Super_Subject_ID)
  g_test_data <- total_data %>% filter(Super_Subject_ID %in% hold_out_test$Super_Subject_ID)
  g_full_train_data <- total_data %>% filter(Super_Subject_ID %in% full_train_data$Super_Subject_ID)
  
  
  # Extract dataset name only by removing .csv from the dataset file name
  dataset_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
  
  # Partitioning hold_out_test data further into Testset 1a, Testset 1b or  Testset 2a, Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072" | dataset_name == "Gene_Expression_Dataset_2_GSE68310") {
    set.seed(data_partition_seed)
    # Dividing test data set further into holdout test set (25%) and validation set (25%) using createDataPartition function of caret package
    index_test1 <- createDataPartition(y = hold_out_test$True_Class_Label,
                                       p = 0.50, list = FALSE)
    hold_out_test_a <- hold_out_test[index_test1, ]
    hold_out_test_b <- hold_out_test[-index_test1, ]
    
    # Holdout_Testset_a and Holdout_Testet_b all time points
    holdout_test_a <- total_data %>% filter(Super_Subject_ID %in% hold_out_test_a$Super_Subject_ID)
    holdout_test_b <- total_data %>% filter(Super_Subject_ID %in% hold_out_test_b$Super_Subject_ID)
    
    # Return training, validation and holdout testsets for Dataset 1 and Dataset 2
    return(list(train_data = g_train_data, full_train_data = g_full_train_data, valid_data = g_valid_data, holdout_test = g_test_data, holdout_test_a = holdout_test_a, holdout_test_b = holdout_test_b))
    
  } else {
    # Return training, validation and holdout testsets for Dataset 3 and Dataset 4 (Can't partition these datasets into holdout testset a and b as they are small in size)
    return(list(train_data = g_train_data, full_train_data = g_full_train_data, valid_data = g_valid_data, holdout_test = g_test_data))
  }
  
}

# Train KNN model
# If k_vals is specified, train using this, which might be a single value or a set of values for a grid search. Otherwise, use default settings for a grid search.
train_knn_model <- function(train_data, Labels_train_data, ml_model_seed, k_vals=seq(1,dim(train_data)[1],2)) {
  set.seed(ml_model_seed)
  metric <- "Accuracy"
  
  grid <- expand.grid(k = k_vals)
  trained_model <- train(x= train_data,
                         y = Labels_train_data,
                         method = "knn",
                         metric = metric,
                         tuneGrid = grid)
  
  return(trained_model)
}


# Train Random Forest model 
train_rf_model <- function(train_data, Labels_train_data, ml_model_seed, mtry = seq(10,floor(sqrt(ncol(train_data))),10), ntree_vals = 100) {
  
  # defining evaluation metric
  metric <- "Accuracy"
  # ntree: parameter that allows number of trees to grow
  # The mtry parameter setting: Number of variables selected as candidates at each split.
  
  # Passing parameter into tunegrid
  grid <- expand.grid(.mtry=mtry)
  set.seed(ml_model_seed)
  trained_model <- train(x= train_data,
                         y = Labels_train_data,
                         method = "rf",
                         metric = metric,
                         tuneGrid = grid,
                         ntree = ntree_vals)
  return(trained_model)
}


# Train LSVM model 
train_LSVM_model <- function(train_data, Labels_train_data, ml_model_seed, 
                             c_vals=c(2^-5, 2^-3, 2^-1, 1, 2^1, 2^3, 2^5,
                                      2^7, 2^9, 2^11, 2^13, 2^15)) {
  set.seed(ml_model_seed)
  
  # Assigning values to the parameter C
  grid <- expand.grid(C = c_vals)
  
  # training LSVM classifier
  trained_model <- train(x= train_data, # Training Data
                         y = Labels_train_data,  # Class labels of training data
                         method = "svmLinear", # Train using LSVM
                         tuneGrid = grid) # Passing training control parameters
  # Print trained model
  return(trained_model)
}


# Train RBF SVM model 
train_RBF_SVM_model <- function(train_data, Labels_train_data, ml_model_seed,
                                c_vals = c(2^-5, 2^-3, 2^-1, 1, 2^1, 2^3, 2^5,
                                           2^7, 2^9, 2^11, 2^13, 2^15),
                                sigma_vals = c(2^-15, 2^-13, 2^-11, 2^-9, 2^-7, 2^-5, 2^-3, 
                                               2^-1, 1, 2^1, 2^3)) {
  set.seed(ml_model_seed)
  
  # Assigning values to the parameter sigma and C
  gridRBF <- expand.grid(sigma = sigma_vals,
                         C = c_vals)
  
  # training LSVM classifier
  trained_model <- train(x= train_data, # Training Data
                         y = Labels_train_data,  # Class labels of training data
                         method = "svmRadial", # Train using RBF SVM
                         tuneGrid = gridRBF) # Passing training control parameters
  # Print trained model
  return(trained_model)
}


# Train XGBoost model 
train_XGBoost_model <- function(train_data, Labels_train_data, ml_model_seed, max_depth_vals = c(2:6), eta_vals = seq(0.1,1,0.1), nrounds_vals = seq(10,100,10)) {
  set.seed(ml_model_seed)
  
  # Assigning values to the parameters
  tune_grid <- expand.grid(
    nrounds = nrounds_vals,
    eta = eta_vals,
    max_depth = max_depth_vals,
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  
  # training XGBoost classifier
  trained_model <- caret::train(
    x = train_data,
    y = factor(Labels_train_data),
    tuneGrid = tune_grid,
    method = "xgbTree",
    metric = "Accuracy",
    verbose = FALSE,
    verbosity = 0
  )
  
  
  # Print trained model
  return(trained_model)
}



# Predict function 
predict_test <- function(final_trained_model, holdout_test) {
  
  # predict subject's health at target time point 
  holdout_test <- holdout_test %>% filter(Time > 0)
  
  predictions <- predict(final_trained_model, newdata = as.matrix(holdout_test[,-c(1:6)]))
  
  # return predictions
  return(predictions)
}

# Write confusion matrix results to TXT
write_confusion_to_txt <- function(predictions, holdout_test, model_name, dataset_name, output_path, KB_name = "____") {
  
  # predict subject's health at target time point 
  holdout_test <- holdout_test %>% filter(Time > 0)
  
  actual_labels <- as.factor(holdout_test$True_Class_Label)
  
  # Ensure that both predictions and actual_labels are factors
  predictions <- factor(predictions)
  actual_labels <- factor(actual_labels)
  
  dataset_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
  KB_name <- substr(KB_name, 1, nchar(KB_name) - 4)
  # Create the filename
  filename <- paste0(output_path, model_name, "_", dataset_name,"_", KB_name, "_confusion_matrix.txt")
  
  # Open a connection for writing
  con <- file(filename, "w")
  print(paste0("Writing results to ", output_path))
  
  if(length(levels(predictions)) == 1 && length(levels(actual_labels)) == 1) {
    # Calculate accuracy for a single class
    accuracy <- sum(predictions == actual_labels) / length(actual_labels)
    write("Only one class present in both predictions and actual labels. Calculated accuracy for single class.", con)
    write(print(paste0("Accuracy: ", accuracy*100,"%")), con)
    
    print("Only one class present in both predictions and actual labels.")
    print(paste0("Calculated accuracy for single class. Accuracy: ", accuracy*100,"%"))
  } else {
    # Ensure both have the same levels
    levels_union <- union(levels(predictions), levels(actual_labels))
    predictions <- factor(predictions, levels = levels_union)
    actual_labels <- factor(actual_labels, levels = levels_union)
    
    # Compute the confusion matrix
    conf_matrix_result <- confusionMatrix(predictions, actual_labels)
    
    # Extract the matrix and statistics
    conf_matrix_table <- conf_matrix_result$table
    overall_stats <- conf_matrix_result$overall
    class_stats <- conf_matrix_result$byClass
    
    # Print and write confusion matrix
    print(paste0("Print confusion matrix for hold-out test set for ", model_name, " ", "for ", dataset_name, ":"))
    print(conf_matrix_result)
    
    write("Confusion Matrix:", con)
    write.table(conf_matrix_table, con, sep = "\t")
    
    write("\nOverall Statistics:", con)
    write.table(as.data.frame(overall_stats), con, sep = "\t", row.names = TRUE)
    
    write("\nClass Statistics:", con)
    write.table(as.data.frame(t(class_stats)), con, sep = "\t", row.names = TRUE)
  }
  
  # Close the connection
  close(con)
}


##################################################################################################
####################################### Main Program  ############################################
##################################################################################################

# Read datasets
dataset_names <- c("Gene_Expression_Dataset_4_GSE61754.csv", "Gene_Expression_Dataset_3_GSE90732.csv", "Gene_Expression_Dataset_2_GSE68310.csv", "Gene_Expression_Dataset_1_GSE73072.csv")


for(dataset_name in dataset_names) {
  Gene_Exp_Data <- read_gene_expression_data(dataset_name, input_path)
  
  # Split data
  splits <- split_data(Gene_Exp_Data, data_partition_seed, dataset_name)
  
  ##################################################################################################
  ########################################### KNN  #################################################
  ##################################################################################################
  
  # Model building using training data
  print("########### Starting KNN learning using training data ###########")
  
  # Model building using training data
  trained_knn_model <- train_knn_model(as.matrix(splits$train_data[,-c(1:6)]), # supply train data after removing labels and other demographic info columns. Only supply gene expression data
                                       as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print KNN training results
  print("KNN training results:")
  print(trained_knn_model)
  
  # Performing validation and hyper parameter selection using validation data
  validation_knn_model <- train_knn_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                          as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  
  # Print KNN validation results
  #print(paste0("KNN validation results: ", validation_knn_model))
  cat("\n")
  print("########### Starting validation process ###########")
  cat("\n")
  print("KNN validation results:")
  print(validation_knn_model)
  
  
  # Selecting final model parameters
  final_k <- validation_knn_model$finalModel$tuneValue[1]
  
  # Print final parameter values
  print(paste0("Final value of k: ", final_k))
  
  # Final model building using KNN
  final_knn_trained_model <- train_knn_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                             as.factor(splits$full_train_data$True_Class_Label),
                                             ml_model_seed, final_k)
  
  # Test KNN
  cat("\n")
  print("########### Starting prediction for holdout testset ###########")
  cat("\n")
  knn_predictions <- predict_test(final_knn_trained_model, splits$holdout_test)
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(knn_predictions, splits$holdout_test, "KNN", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test KNN
    knn_predictions <- predict_test(final_knn_trained_model, splits$holdout_test_a)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file for holdout_testset_a
    write_confusion_to_txt(knn_predictions, splits$holdout_test_a, "KNN", result_filename, output_path)
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test KNN
    knn_predictions <- predict_test(final_knn_trained_model, splits$holdout_test_b)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file for holdout_testset_b
    write_confusion_to_txt(knn_predictions, splits$holdout_test_b, "KNN", result_filename, output_path)
  }
  
  ##################################################################################################
  ######################################## End KNN  ################################################
  ##################################################################################################
  
  
  ##################################################################################################
  ##########################################  RF  ##################################################
  ##################################################################################################
  
  # Code for RF model building, validation and evaluation
  # Model building using training data
  print("########### Starting RF learning using training data ###########")
  trained_rf_model <- train_rf_model(as.matrix(splits$train_data[,-c(1:6)]),
                                     as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print RF training results
  print("RF training results:")
  print(trained_rf_model)
  
  # Performing validation and hyper parameter selection using validation data
  validation_rf_model <- train_rf_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                        as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  # Print RF validation results
  print("########### RF validation results ###########")
  print(validation_rf_model)
  
  # Selecting final model parameters
  final_mtry <- validation_rf_model$finalModel$mtry
  final_n_tree <- validation_rf_model$finalModel$ntree
  
  # Print final parameter values
  print("Final value of mtry (na) is:")
  print(final_mtry)
  print("Final value of n_tree (nt) is:")
  print(final_n_tree)
  
  
  # Final model building using RF
  final_rf_trained_model <- train_rf_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                           as.factor(splits$full_train_data$True_Class_Label),
                                           ml_model_seed, final_mtry, final_n_tree)
  
  # Test RF model
  print("########### Starting prediction for holdout testset ###########")
  rf_predictions <- predict_test(final_rf_trained_model, splits$holdout_test)
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(rf_predictions, splits$holdout_test, "RF", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test RF
    rf_predictions <- predict_test(final_rf_trained_model, splits$holdout_test_a)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file for holdout_testset_a
    write_confusion_to_txt(rf_predictions, splits$holdout_test_a, "RF", result_filename, output_path)
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test RF
    rf_predictions <- predict_test(final_rf_trained_model, splits$holdout_test_b)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file for holdout_testset_b
    write_confusion_to_txt(rf_predictions, splits$holdout_test_b, "RF", result_filename, output_path)
  }
  
  
  ##################################################################################################
  ########################################  End RF  ################################################
  ##################################################################################################
  
  
  ##################################################################################################
  ####################################### Linear SVM ###############################################
  ##################################################################################################
  
  # Code for LSVM model building, validation and evaluation
  # Model building using training data
  print("###########  Starting Linear SVM learning using training data ########### ")
  trained_LSVM_model <- train_LSVM_model(as.matrix(splits$train_data[,-c(1:6)]),
                                         as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print LSVM training results
  print("LSVM training results:")
  print(trained_LSVM_model)
  
  # Performing validation and hyper parameter selection using validation data
  validation_LSVM_model <- train_LSVM_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                            as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  # Print LSVM validation results
  print("########### LSVM validation results ########### ")
  print(validation_LSVM_model)
  
  # Selecting final model parameters
  final_C <- validation_LSVM_model$finalModel@param$C
  
  
  # Print final parameter values
  print("Final value of parameter C of LSVM:")
  print(final_C)
  
  # Final model building using LSVM
  final_LSVM_trained_model <- train_LSVM_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                               as.factor(splits$full_train_data$True_Class_Label),
                                               ml_model_seed, final_C)
  
  # Test LSVM model
  print("########### Starting prediction for holdout testset ###########")
  LSVM_predictions <- predict_test(final_LSVM_trained_model, splits$holdout_test)
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(LSVM_predictions, splits$holdout_test, "LSVM", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test LSVM
    LSVM_predictions <- predict_test(final_LSVM_trained_model, splits$holdout_test_a)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file for holdout_testset_a
    write_confusion_to_txt(LSVM_predictions, splits$holdout_test_a, "LSVM", result_filename, output_path)
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test LSVM
    LSVM_predictions <- predict_test(final_LSVM_trained_model, splits$holdout_test_b)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file for holdout_testset_b
    write_confusion_to_txt(LSVM_predictions, splits$holdout_test_b, "LSVM", result_filename, output_path)
  }
  
  
  ##################################################################################################
  ######################################  End Linear SVM  ##########################################
  ##################################################################################################
  
  
  ##################################################################################################
  ######################################### RBF SVM ################################################
  ##################################################################################################
  
  # Code for RBF_SVM model building, validation and evaluation
  # Model building using training data
  print("###########  Starting RBF SVM learning using training data ########### ")
  trained_RBF_SVM_model <- train_RBF_SVM_model(as.matrix(splits$train_data[,-c(1:6)]),
                                               as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print RBF_SVM training results
  print("RBF SVM training results:")
  print(trained_RBF_SVM_model)
  
  # Performing validation and hyper parameter selection using validation data
  validation_RBF_SVM_model <- train_RBF_SVM_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                                  as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  # Print RBF_SVM validation results
  print("###########  RBF SVM validation results ########### ")
  print(validation_RBF_SVM_model)
  
  # Selecting final model parameters
  final_sigma <- validation_RBF_SVM_model$bestTune$sigma
  final_C = validation_RBF_SVM_model$bestTune$C
  
  
  # Print final parameter values
  print("Final value of parameter C of RBF SVM:")
  print(final_C)
  print("Final value of parameter sigma of RBF SVM:")
  print(final_sigma)
  
  # Final model building using SVM
  final_RBF_SVM_trained_model <- train_RBF_SVM_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                                     as.factor(splits$full_train_data$True_Class_Label),
                                                     ml_model_seed, final_C, final_sigma)
  
  # Test RBF_SVM model
  print("########### Starting prediction for holdout testset ###########")
  SVM_predictions <- predict_test(final_RBF_SVM_trained_model, splits$holdout_test)
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(SVM_predictions, splits$holdout_test, "RBF_SVM", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test RBF_SVM
    SVM_predictions <- predict_test(final_RBF_SVM_trained_model, splits$holdout_test_a)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file for holdout_testset_a
    write_confusion_to_txt(SVM_predictions, splits$holdout_test_a, "RBF_SVM", result_filename, output_path)
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    # Test RBF_SVM
    SVM_predictions <- predict_test(final_RBF_SVM_trained_model, splits$holdout_test_b)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file for holdout_testset_b
    write_confusion_to_txt(SVM_predictions, splits$holdout_test_b, "RBF_SVM", result_filename, output_path)
  }
  
  
  ##################################################################################################
  ########################################  End RBF SVM  ###########################################
  ##################################################################################################
  
  
  ##################################################################################################
  ######################################### XGBoost ################################################
  ##################################################################################################
  
  # Code for XGBoost model building, validation and evaluation
  # Model building using training data
  print("###########  Starting XGBoost learning using training data ########### ")
  trained_XGBoost_model <- train_XGBoost_model(as.matrix(splits$train_data[,-c(1:6)]),
                                               as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print XGBoost training results
  print("XGBoost training results:")
  print(trained_XGBoost_model)
  
  # Performing validation and hyper parameter selection using validation data
  validation_XGBoost_model <- train_XGBoost_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                                  as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  # Print XGBoost validation results
  print("########### XGBoost validation results ###########")
  print(validation_XGBoost_model)
  
  # Selecting final model parameters
  final_max_depth <- validation_XGBoost_model$bestTune$max_depth
  final_eta <- validation_XGBoost_model$bestTune$eta
  final_nrounds <- validation_XGBoost_model$bestTune$nrounds
  
  # Print final parameter values
  print("Final value of parameter max depth of XGBoost:")
  print(final_max_depth)
  print("Final value of parameter eta of XGBoost:")
  print(final_eta)
  print("Final value of parameter rounds of XGBoost:")
  print(final_nrounds)
  
  
  # Final model building using XGBoost
  final_XGBoost_trained_model <- train_XGBoost_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                                     as.factor(splits$full_train_data$True_Class_Label),
                                                     ml_model_seed, final_max_depth, final_eta, final_nrounds)
  
  # Test XGBoost model
  print("########### Starting prediction for holdout testset ###########")
  XGBoost_predictions <- predict_test(final_XGBoost_trained_model, splits$holdout_test)
  
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(as.factor(XGBoost_predictions), splits$holdout_test, "XGBoost", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    
    # Test XGBoost model
    XGBoost_predictions <- predict_test(final_XGBoost_trained_model, splits$holdout_test_a)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file
    write_confusion_to_txt(as.factor(XGBoost_predictions), splits$holdout_test_a, "XGBoost", result_filename, output_path)
    
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    
    # Test XGBoost model
    XGBoost_predictions <- predict_test(final_XGBoost_trained_model, splits$holdout_test_b)
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file
    write_confusion_to_txt(as.factor(XGBoost_predictions), splits$holdout_test_b, "XGBoost", result_filename, output_path)
    
  } # End if condition
  
  ##################################################################################################
  ########################################  End XGBoost  ###########################################
  ##################################################################################################
  
  
}
