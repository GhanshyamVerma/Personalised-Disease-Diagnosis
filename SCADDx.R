# Required Programming Language
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

# Function to read Knowledge Bases
read_KBs <- function(filename,input_path) {
  path_filename <- paste0(input_path, filename)
  print(paste0("Reading KB called ", path_filename))
  return(read.csv(file = path_filename, header = TRUE, sep=","))
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



# Train SCADDx model 
SCADDx_model <- function(KB, train_data, Labels_train_data, ml_model_seed, dataset_name, KB_name, output_path, model_name, which_set, P_start = 25 , Q_start = 25) {
  set.seed(ml_model_seed)
  
  # enter the number of top "m" most likely diseases you want to see with computed probabilities for each patient 
  m = 5
  
  # enter values of MDEGs and LDEGs you want to consider or enter a range
  # if you don't want to perform grid search and want to perform the experiment on single P and Q value 
  # then assign same value to variable P_start and P_end and assign 0 value to P_step = 0. Similarly same value (value of Q) to variable Q_start and Q_end and value 0 to Q_step as shown below.
  # P_start <- 25
  # P_end <- 25
  # P_step <- 0
  # Q_start <- 25
  # Q_end <- 25
  # Q_step <- 0
  
  
  if(which_set == "Holdout_test" | which_set == "Final_model"){
    ###################################    NO Grid search   ##############################################################
    
    P_start <- P_start
    P_end <- P_start
    P_step <- 0
    Q_start <- Q_start
    Q_end <- Q_start
    Q_step <- 0
    
  }else{
    ###################################    Grid search   ##############################################################
    
    # Grid search: alternatively to perform grid search, enter appropriate start, end, and step size for P and Q
    P_start <- 25
    P_end <- 50      # replace P_end <- 50 with 300 to reproduce the results shown in the paper. Keep it low for quick test.  
    P_step <- 25
    Q_start <- 25
    Q_end <- 50      # replace Q_end <- 50 with 300 to reproduce the results shown in the paper. Keep it low for quick test.
    Q_step <- 25
    
  }
  
  
  #######################################################################################################################################
  
  print(paste0("##################################################### ",which_set," #############################################################"))
  
  #######################################################################################################################################
  
  
  # extract all unique diseases from KB to assign computed disease weight
  Data_Unique_Disease <- KB[!duplicated(KB[,c('disease_id')]), c('disease_name','disease_id')]
  
  # reorder the columns
  Data_Unique_Disease <- Data_Unique_Disease[ , c('disease_id','disease_name')]
  
  # add a new column named Disease_Weight
  Data_Unique_Disease <- Data_Unique_Disease %>% mutate(Disease_Weight = 0)
  
  # temporary variable of Data_Unique_Disease
  Data_Unique_Disease_initial_weights <- Data_Unique_Disease
  
  # gene expression values start index in gene expression data
  s_index <- 7
  
  # computing length for grid search
  PS <- length(seq(P_start,P_end,P_step))
  QS <- length(seq(Q_start,Q_end,Q_step))
  
  # creating data frame to compute and store accuracy at different values of P and Q and, and different value of top n diseases
  Accuracy_matrix <- data.frame("P" = 1:(PS*QS), "Q" = 1:(PS*QS), "Acc_Top_1_Dis" = 1:(PS*QS), "Acc_Top_2_Dis" = 1:(PS*QS), "Acc_Top_3_Dis" = 1:(PS*QS), "Acc_Top_4_Dis" = 1:(PS*QS), "Acc_Top_5_Dis" = 1:(PS*QS), "Acc_Top_10_Dis" = 1:(PS*QS))
  
  Accuracy_matrix <- Accuracy_matrix %>% mutate(Avg_Acc = 0)
  Accuracy_matrix <- Accuracy_matrix %>% mutate(Time = 0)
  # initializing the incrementer for Accuracy_matrix
  Acc_index <- 1
  
  # Extract dataset name only by removing .csv from the dataset file name
  dataset_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
  
  # Extract KB name only by removing .csv from the KB file name
  KB_name <- substr(KB_name, 1, nchar(KB_name) - 4)
  
  # Start the clock!
  start_loop_time <- proc.time()
  
  #### Code for assigning weights to the diseases in KB based on changes observed in genes
  
  for(p in seq(P_start,P_end,P_step)){  # loop of p is for top genes / P MDEGs
    
    for(q in seq(Q_start,Q_end,Q_step)){  # loop of q is for bottom genes / Q LDEGs
      
      # Start the clock!
      start_p_time <- proc.time()
      
      # total number of subjects
      s <- dim(train_data)[1]/2
      
      # making super_subject_id serial wise
      train_data$Super_Subject_ID <- rep(1:s, each = 2)
      
      # creating data frame to store values of predicted class labels
      
      predicted_info <- train_data[ , c(1:s_index-1)]
      predicted_info <- predicted_info %>% mutate(predicted_label_top_1 = 1:(2*s))
      predicted_info <- predicted_info %>% mutate(predicted_label_top_2 = 1:(2*s))
      predicted_info <- predicted_info %>% mutate(predicted_label_top_3 = 1:(2*s))
      predicted_info <- predicted_info %>% mutate(predicted_label_top_4 = 1:(2*s))
      predicted_info <- predicted_info %>% mutate(predicted_label_top_5 = 1:(2*s))
      predicted_info <- predicted_info %>% mutate(predicted_label_top_10 = 1:(2*s))
      
      Gene_Data_All_ti_prediction <- train_data[!train_data$Time == 0, c(1:(s_index-1))]
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_1 = 1:s)
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_2 = 1:s)
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_3 = 1:s)
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_4 = 1:s)
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_5 = 1:s)
      Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_10 = 1:s)
      
      All_Sub_temp_prediction <- data.frame("Top_1"=1:s, "Top_2"=1:s,"Top_3"=1:s, "Top_4"=1:s,"Top_5"=1:s, "Top_10"=1:s)
      
      # loop of l for number of subjects
      
      for(l in 1:s){ 
        
        # extracting data of lth subject 
        # making super_subject_id serial wise
        # train_data$Super_Subject_ID <- rep(1:s, each = 2)
        Gene_expression_data_sub_l <- train_data %>% filter(Super_Subject_ID == l)
        Gene_expression_data_sub_l <- Gene_expression_data_sub_l[ , -c(1:(s_index-1))]
        
        cat("\n")
        print(paste0("########################### ",which_set," code: New subject's computation start from here#############################"))
        
        cat("\n")
        print(paste0(which_set," Data: Subject/Patient id is:"))
        print(l)
        cat("\n")
        
        # for each subject, again initialize the Data_Unique_Disease variable dataframe
        Data_Unique_Disease <- Data_Unique_Disease_initial_weights
        
        # compute changes in gene expression values (Target sample (TD) - Reference sample (T1)
        Gene_Transition_Matrix <- Gene_expression_data_sub_l[2, ] - Gene_expression_data_sub_l[1, ]
        
        # extracting top P MDEGs
        Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix[ , order(-abs(Gene_Transition_Matrix[ , ]))]
        Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix_top_p_Genes[ , c(1:p)]
        
        print("Top 5 Most Differencially Expressed Genes (MDEGs):")
        print(Gene_Transition_Matrix_top_p_Genes[ , 1:5])
        cat("\n")
        
        for(j in 1:p){ # loop of j for number of genes
          
          # extract all diseases' ids from KB that are associated with top P genes/ MDEGs of the subject
          Disease_IDs <- KB[KB$gene_symbol == names(Gene_Transition_Matrix_top_p_Genes)[j], "disease_id" ]
          
          if(length(Disease_IDs) == 0){
          }else{
            for(k in 1:length(Disease_IDs)){ # loop for every disease id
              
              # computing the weight/score for each disease for lth subject
              Data_Unique_Disease[Data_Unique_Disease$disease_id ==  Disease_IDs[k], "Disease_Weight"] <-   Data_Unique_Disease[Data_Unique_Disease$disease_id ==  Disease_IDs[k], "Disease_Weight"] + abs(Gene_Transition_Matrix_top_p_Genes[,j])
              
            }  # end for loop k
            
          } # end else
          
          
        } # end for loop j
        
        #### Code for down-weighting the diseases based on the associated Q LDEGs to them
        
        # extracting Q LDEGs
        Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix[ , order(abs(Gene_Transition_Matrix[ , ]))]
        Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix_bottom_q_Genes[ , c(1:q)]
        
        print("5 Least Differencially Expressed Genes (LDEGs):")
        print(Gene_Transition_Matrix_bottom_q_Genes[ , 1:5])
        cat("\n")
        
        for(j in 1:q){ # loop of j for number of bottom genes
          
          # extract all diseases' ids from KB that are associated with bottom Q genes / LDEGs of the subject
          Disease_IDs <- KB[KB$gene_symbol == names(Gene_Transition_Matrix_bottom_q_Genes)[j], "disease_id" ]
          
          if(length(Disease_IDs) == 0){
          }else{
            for(k in 1:length(Disease_IDs)){ # loop for every disease id
              
              # down-weighting the diseases based on the associated Q LDEGs to them
              Data_Unique_Disease[Data_Unique_Disease$disease_id ==  Disease_IDs[k], "Disease_Weight"] <-  Data_Unique_Disease[Data_Unique_Disease$disease_id ==  Disease_IDs[k], "Disease_Weight"] - (abs(Gene_Transition_Matrix_top_p_Genes[,1]) - abs(Gene_Transition_Matrix_bottom_q_Genes[,j]))
              
            }  # end for loop k
            
          } # end else
          
          
        } # end for loop j
        
        
        # create file name to write data into csv file
        #file_name <- paste("Disease_Weight_Train_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
        
        # write data into csv file
        #write.csv(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ], file = file_name, row.names = FALSE)
        
        print("Value of P (number of MDEGs) is :")
        print(p)
        cat("\n")
        
        print("Value of Q is :")
        print(q)
        cat("\n")
        
        print("This subject at this time point has following True Class Label:")
        print(train_data[ train_data$Super_Subject_ID == l , ]$True_Class_Label[2])
        cat("\n")
        
        for(i in 1:6){ # loop for how many top disease you want to look for Acc calc (current loop is for top 1 to 5 and top 10 predicted diseases )
          
          if(i<6){
            Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:i, "disease_name"]
          }else{
            Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, "disease_name"]
          }
          
          if(any(Top_Disease_Names == "Respiratory Viral Infection") || any(Top_Disease_Names == "Influenza") || any(Top_Disease_Names == "Respiratory Tract Diseases") || any(Top_Disease_Names == "Respiratory Syncytial Virus Infections") || any(Top_Disease_Names == "Respiratory Tract Infections")){
            predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "RVI"
            Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "RVI"
            All_Sub_temp_prediction[l,i] <- "RVI"
          }else{
            predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "Not RVI"
            Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "Not RVI"
            All_Sub_temp_prediction[l,i] <- "Not RVI"
          }
          
          if(i<6){
            print(paste("Predicted label using top ", i, "disease is:"))
            print(All_Sub_temp_prediction[l,i])
            cat("\n")
            
          }else{
            print("Predicted label using top 10 disease is:")
            print(All_Sub_temp_prediction[l,i])
            cat("\n")
            
          }
          
        }
        
        # creating dataframe to store computed disease probabilites
        Data_Unique_Disease_with_Probability <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:m, ]
        
        # function to compute probability of predicted disease using generalized softmax
        generalized_softmax <- function(Dis_weight){
          S_Prob <- Dis_weight
          e_Dis_weight <- exp(Dis_weight)
          for(i in 1:length(Dis_weight)){
            S_Prob[i] <- (exp(Dis_weight[i])/sum(e_Dis_weight))*100
          }
          return(S_Prob)
        }
        
        # calling probability function
        Disease_Prob <- generalized_softmax(Data_Unique_Disease_with_Probability$Disease_Weight)
        Data_Unique_Disease_with_Probability <- Data_Unique_Disease_with_Probability %>% mutate(Disease_Probability = Disease_Prob)
        
        # extracting other entities from the KB for the top predicted disease
        disease_class_n_others <- KB[1:m, c('disease_type', 'disease_class_name', 'disease_semantic_type','gene_symbol', 'protein_class', 'uniprot_id')]
        disease_class_n_others <- disease_class_n_others %>% mutate(change_in_gene_expr = 0)
        
        for(d in 1:m){
          disease_class_n_others$disease_type[d] <- KB[KB$disease_id == Data_Unique_Disease_with_Probability$disease_id[d], c('disease_type')][1]
          disease_class_n_others$disease_class_name[d] <- KB[KB$disease_id == Data_Unique_Disease_with_Probability$disease_id[d], c('disease_class_name')][1] 
          disease_class_n_others$disease_semantic_type[d] <- KB[KB$disease_id == Data_Unique_Disease_with_Probability$disease_id[d], c('disease_semantic_type')][1] 
          disease_class_n_others$gene_symbol[d] <- colnames(Gene_Transition_Matrix_top_p_Genes)[d]
          disease_class_n_others$change_in_gene_expr[d] <- Gene_Transition_Matrix_top_p_Genes[,d]
          disease_class_n_others$protein_class[d] <- KB[KB$gene_symbol == colnames(Gene_Transition_Matrix_top_p_Genes)[d] & KB$disease_id == Data_Unique_Disease_with_Probability$disease_id[1], c('protein_class')][1] 
          disease_class_n_others$uniprot_id[d] <- KB[KB$gene_symbol == colnames(Gene_Transition_Matrix_top_p_Genes)[d] & KB$disease_id == Data_Unique_Disease_with_Probability$disease_id[1], c('uniprot_id')][1] 
          
        }
        
        Data_Unique_Disease_with_Probability <- cbind(Data_Unique_Disease_with_Probability,disease_class_n_others)
        
        
        print("Top 5 most likely diseases predicted for this subject are:")
        print(Data_Unique_Disease_with_Probability)
        cat("\n")
        
        cat("\n")
        print("################################## This subject's predictions ends here #####################################")
        cat("\n")
      } # end for loop 
      
      
      if(which_set == "Holdout_test"){
        # Create the file name
        filename <- paste0(output_path, model_name, "_", dataset_name,"_",KB_name,"_predicted_info_",which_set,"_Total_Sub",l,"_p_",p,"_q_",q,".csv")
        write.csv(Gene_Data_All_ti_prediction, file = filename, row.names = FALSE)
      }
      
      
      cat("\n")
      print("Accuracy calculated using all subjects considering different values of top n diseases:")
      cat("\n")
      
      ## computing accuracy using all subjects of test data 
      # computing accuracy using confusion matrix if both positive and negative subjects are there in the data
      if(any(Gene_Data_All_ti_prediction$True_Class_Label == "Not RVI")){
        for(i in 1:6){
          confusion_mat <- confusionMatrix( as.factor(All_Sub_temp_prediction[,i]), as.factor(Gene_Data_All_ti_prediction$True_Class_Label), positive = "RVI")
          
          if(i<6){
            cat("\n")
            print(paste("Accuracy using top", i, "disease is:"))
            print(confusion_mat)
            cat("\n")
            
          }else{
            cat("\n")
            print("Accuracy using top 10 diseases is:")
            print(confusion_mat)
            cat("\n")
            
          }
          
          Accuracy_matrix[Acc_index, (2+i)] <- confusion_mat$overall[1]
        }
      }else{ #computing accuracy using standard method (hit rate) if only positive subjects are there in the data
        for(i in 1:6){
          hit <- 0
          for(k in 1:dim(All_Sub_temp_prediction)[1]){
            if(Gene_Data_All_ti_prediction[k,"True_Class_Label"] == All_Sub_temp_prediction[k,i]){
              hit <- hit + 1
            }
          }
          Acc <- hit/dim(All_Sub_temp_prediction)[1]
          cat("\n")
          print(paste("Accuracy using top", i, "disease is:"))
          print(Acc)
          if(i<6){
            cat("\n")
            print(paste("Accuracy using top", i, "disease is:"))
            print(Acc)
            cat("\n")
            
          }else{
            cat("\n")
            print("Accuracy using top 10 diseases is:")
            print(Acc)
            cat("\n")
            
          }
          Accuracy_matrix[Acc_index, (2+i)] <- Acc
          
        }
        
      }
      
      # assign current value of P and Q 
      Accuracy_matrix$P[Acc_index] <- p
      Accuracy_matrix$Q[Acc_index] <- q
      
      
      
      # Stop the clock
      each_itteration_time <- proc.time() - start_p_time
      
      print(each_itteration_time)
      # adding time to the Accuracy_matrix
      Accuracy_matrix$Time[Acc_index] <- each_itteration_time[3]
      
      # adding average acc to the Accuracy_matrix
      Accuracy_matrix$Avg_Acc <- rowMeans(Accuracy_matrix[,c(3:8)])
      
      cat("\n")
      print("Accuracy calculated using all subjects on different values of top n diseases and varing P and Q:")
      cat("\n")
      print("Accuracy Matrix:")
      print(Accuracy_matrix[1:Acc_index, ])
      cat("\n")
      
      
      # incrementer 
      Acc_index <- Acc_index +1
      
      
      
      
      print("################################## Current itteration of the loop for Q ends here ##########################################")
      cat("\n")
      
    } # ending loop q 
    
    print("################################## Current itteration of the loop for P ends here ##########################################")
    cat("\n")
    
  } # ending loop p
  
  # Stop the clock
  total_loop_time <- proc.time() - start_loop_time
  
  print(total_loop_time)
  
  if(which_set == "Holdout_test"){
    # Create the filename
    filename <- paste0(output_path, model_name, "_", dataset_name,"_",KB_name,"_Accuracy_Matrix_",which_set,"_Total_Sub",l,"_p_",p,"_q_",q,".csv")
    write.csv(Accuracy_matrix, file = filename, row.names = FALSE)
    
  }
  
  
  
  return(list(Accuracy_matrix = Accuracy_matrix, Predictions = Gene_Data_All_ti_prediction))
  
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
  
  # Read Knowledge Bases
  
  KB_names <- c("DisGeNet_Knowledge_Base.csv", "CTD_Knowledge_Base.csv")
  for(KB_name in KB_names) {
    KB <- read_KBs(KB_name, input_path)
    
    
    print(paste0("Dataset name: ", dataset_name))
    print(paste0("KB name: ", KB_name))
    
    ##################################################################################################
    ######################################## SCADDx  #################################################
    ##################################################################################################
    
    # Model building using training data
    print("########### Starting SCADDx learning using training data ###########")
    
    trained_SCADDx_model <- SCADDx_model(KB, splits$train_data, splits$train_data$True_Class_Label, ml_model_seed, dataset_name, KB_name, output_path, "SCADDx", "Train")
    
    
    # Print SCADDx training results
    print("SCADDx training results:")
    print(trained_SCADDx_model$Accuracy_matrix)
    
    # Performing validation and hyper parameter selection using validation data
    validation_SCADDx_model <- SCADDx_model(KB, splits$valid_data, splits$valid_data$True_Class_Label, ml_model_seed, dataset_name, KB_name, output_path, "SCADDx", "Validation")
    
    
    # Print SCADDx validation results
    print("SCADDx validation results:")
    print(validation_SCADDx_model$Accuracy_matrix)
    
    # Selecting best value of P and Q based on the performance on validation set
    P <- validation_SCADDx_model$Accuracy_matrix[validation_SCADDx_model$Accuracy_matrix$Avg_Acc == max(validation_SCADDx_model$Accuracy_matrix$Avg_Acc), "P"][1]
    Q <- validation_SCADDx_model$Accuracy_matrix[validation_SCADDx_model$Accuracy_matrix$Avg_Acc == max(validation_SCADDx_model$Accuracy_matrix$Avg_Acc), "Q"][1]
    
    # Print final parameter values
    print(paste0("Final value of P is: ", P))
    print(paste0("Final value of Q is: ", Q))
    
    # Final model building using SCADDx
    final_SCADDx_trained_model <- SCADDx_model(KB, splits$full_train_data, splits$full_train_data$True_Class_Label, ml_model_seed, dataset_name, KB_name, output_path, "SCADDx", "Final_model", P, Q)
    
    
    # Test SCADDx
    cat("\n")
    print("########### Starting prediction for holdout testset ###########")
    cat("\n")
    
    SCADDx_predictions <- SCADDx_model(KB, splits$holdout_test, splits$holdout_test$True_Class_Label, ml_model_seed, dataset_name, KB_name, output_path, "SCADDx", "Holdout_test", final_SCADDx_trained_model$Accuracy_matrix$P, final_SCADDx_trained_model$Accuracy_matrix$Q)
    
    # Print SCADDx holdout testset results
    print("SCADDx holdout testset results: Accuracy on different values of top n diseases")
    print(SCADDx_predictions$Accuracy_matrix)
    
    # Write results to a text file for holdout testset
    write_confusion_to_txt(SCADDx_predictions$Predictions$predicted_label_top_10, splits$holdout_test, "SCADDx", dataset_name, output_path, KB_name)
    
    
    
    # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
    if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
      
      data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
      # Create the file name
      result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
      
      SCADDx_predictions <- SCADDx_model(KB, splits$holdout_test_a, splits$holdout_test_a$True_Class_Label, ml_model_seed, result_filename, KB_name, output_path, "SCADDx", "Holdout_test", final_SCADDx_trained_model$P, final_SCADDx_trained_model$Accuracy_matrix$Q)
      
      # Print SCADDx holdout testset results
      print(paste0(data_name,"_SCADDx_Holdout_", "Testset_a_results:"))
      print(SCADDx_predictions$Accuracy_matrix)
      
      # Write results to a text file for holdout_testset_a
      write_confusion_to_txt(SCADDx_predictions$Predictions$predicted_label_top_10, splits$holdout_test_a, "SCADDx", result_filename, output_path, KB_name)
      
    }
    
    # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
    if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
      
      data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
      # Create the file name
      result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
      
      SCADDx_predictions <- SCADDx_model(KB, splits$holdout_test_b, splits$holdout_test_b$True_Class_Label, ml_model_seed, result_filename, KB_name, output_path, "SCADDx", "Holdout_test", final_SCADDx_trained_model$Accuracy_matrix$P, final_SCADDx_trained_model$Accuracy_matrix$Q)
      
      # Print SCADDx holdout testset results
      print(paste0(data_name,"_SCADDx_Holdout_", "Testset_b_results:"))
      print(SCADDx_predictions$Accuracy_matrix)
      
      # Write results to a text file for holdout_testset_b
      write_confusion_to_txt(SCADDx_predictions$Predictions$predicted_label_top_10, splits$holdout_test_b, "SCADDx", result_filename, output_path, KB_name)
      
    }
    
    ##################################################################################################
    ####################################### End SCADDx  ##############################################
    ##################################################################################################
    
    
    
  }
  
}
