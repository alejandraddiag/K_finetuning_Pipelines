library(tidyverse)
library(vegan)
library(dplyr)
library(purrr)
library(ggupset)


# Defining K
evenness_K<-function(ASV){
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  shannon_index<-diversity(asv, index="shannon") #Defining shannon index
  S<-specnumber(asv) #Defining Richness
  J<-shannon_index/log(S) #Calculating J
  
  return(J)
}

bray_curtis<-function(ASV){
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  bray<-vegdist(asv, method="bray")
  bray<-as.dist(bray)
  return(bray)
}


K_calculator<-function(ASV){
  #First we make sure everything is numeric
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)
  taxa<-asv%>% #Selecting what are the columns forASV
    select(where(is.numeric))
  N<-nrow(asv) #Number of samples
  
  #Calculating each of the elements
  J<-evenness_K(asv) #Evenness
  BC<-bray_curtis(asv) #Bray Curtis Dissimilarity
  K<-mean(J,na.rm=TRUE)*mean(BC, na.rm=TRUE) #Obtaining K
  return(K)
}

# Function for Weighted K
K_calculator_weighted<-function(ASV, alpha, beta){
  #making sure (as always) everything is numeric!
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  J<-mean(na.omit(evenness_K(asv)))
  BC<-mean(na.omit(bray_curtis(asv)))
  
  K_w<-(J^alpha)*(BC^beta)
  return(K_w)
}

K_weighted_function<-function(ASV, alpha, beta){
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  
  J<-mean(na.omit(evenness_K(asv)))
  BC<-mean(na.omit(bray_curtis(asv)))
  k_values<-c()
  a_i<-c()
  b_i<-c()
  #Iterating through the alpha and beta values
  for (a in alpha){
    for (b in beta){
      K_w<-(J^a)*(BC^b)
      k_values<-append(k_values, K_w)
      #Keeping track of what alpha and beta was use
      a_i<-append(a_i, a)
      b_i<-append(b_i, b)
    }
  }
  return(data.frame(k_values=k_values,
                    a_i,
                    b_i
  ))
}

#Filtering Pipeline I
filtering_ASV<-function(ASV){
  #I define the percentages of the threshold
  min_presence_th<-0.1
  low_abundance_limit<-0.2
  #Make sure all the data is numeric
  asv<-ASV%>%mutate_if(is.integer,as.numeric)%>%
    select(where(is.numeric))
  #Setting the number of samples and ASVs
  num_samples<-nrow(asv)
  num_asv<-ncol(asv)
  #Defining two vectors: one to keep the ASVs and one for bin
  keep<-c()
  bin<-c()
  
  for(i in 1:num_asv){
    seq_count<-asv[,i] #Defining the number of sequences per sample in ASV
    #Filtering paremeters
    non_zero_count<-sum(seq_count!=0) 
    low_value_count<-sum(seq_count<20)
    
    #Filtering the data
    ##First saving all the ASVs that have at least one sample with >20 counts
    if(any(seq_count>20)){
      keep<-c(keep,colnames(asv)[i])
      ##Keep if the ASVs happen to have at least 20 seqs  in the given threshold
    }else if(non_zero_count>=(min_presence_th*num_samples) && low_value_count<(low_abundance_limit*num_samples)){
      keep<-c(keep,colnames(asv)[i])
      ##If the ASV does not have at least 20seq number or be present in at least 2/3 then we discard
    }else{
      bin<-c(bin,colnames(asv)[i])
    }
  }
  return(list(keep = keep, bin = bin))
  
}
# K finetune Pipeline I: Stepwise
k_finetune_cummulative<-function(ASV){
  start.time<-Sys.time()#Getting the start time
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  #Calculating K for raw table
  k_unfiltered<-K_calculator(asv)
  BC_unfiltered<-mean(na.omit(bray_curtis(asv)))
  J_unfiltered<-mean(na.omit(evenness_K(asv)))
  #Running the filtering pipeline
  filtered_asv<-filtering_ASV(asv)
  
  #Applying the filter
  asv0_table<-asv%>%
    select(all_of(filtered_asv$keep))

  #Defining the filtered taxa as a variable
  discarded_asv<-filtered_asv$bin
  sample_id<-rownames(asv0_table)
  
  #Calculating K for filtered table (K_0)
  k_0<-K_calculator(asv0_table)
  BC_0<-mean(na.omit(bray_curtis(asv0_table)))
  J_0<-mean(na.omit(evenness_K(asv0_table)))
  
  #Creating a vector to save the names of the sequentially added ASVs
  added_asv<-c()
  #Vector for K values
  K_values<-c()
  BC_values<-c()
  J_values<-c()
  
  asv_filt_iterative<-asv0_table
  #Iterative sequence for K, BC and J calculation
  for (asv_name in discarded_asv){
    asv_table_i<-cbind(asv_filt_iterative, asv[sample_id,asv_name, drop=FALSE]) ##Appending the discarded taxa
    asv_filt_iterative<-asv_table_i 
    

    ## New K calculation
    k_asv_i<-K_calculator(asv_table_i)
    K_values<-append(K_values, k_asv_i)
    
    ## New BC calculation
    bc_asv_i<-mean(na.omit(bray_curtis(asv_table_i)))
    BC_values<-append(BC_values, bc_asv_i) 
    
    ## New Evennness Calculation
    j_asv_i<-mean(na.omit(evenness_K(asv_table_i)))
    J_values<-append(J_values, j_asv_i)
    
    ##Append the ASV name to the vector of names
    added_asv<-append(added_asv, asv_name)
  }
  
  
  end.time<-Sys.time()
  
  time.take<-round(end.time-start.time,2)
  a<-proc.time()
  #Function output
  return(list(
    final_table=asv_filt_iterative, ##To make sure all filtered taxa have been re-introduced
    asv0_table=asv0_table, ## Filtered table
    values_filtered=data.frame(K_val=k_0,J=J_0, BC=BC_0),
    values_unfiltered=data.frame(K=k_unfiltered,BC=BC_unfiltered, J=J_unfiltered), ## K for the raw table
    df_results=data.frame(K_values=K_values,reintroduced_asv=added_asv, BC_values=BC_values, J_values=J_values), ## List of reintroduced taxa
    time.take=time.take
  ))
  
}

# K combi finetune pipeline: Combinatorial

k_finetune_combi <- function(ASV) {
  start.time <- Sys.time()
  asv <- ASV %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  #K-unfiltered calculation
  k_unfiltered <- K_calculator(asv)
  #BC and J calculations
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_K(asv)))
  
  #Applying the filter
  filtered_asv <- filtering_ASV(ASV)
  
  #Saving the discarded ASVs
  discarded_asv <- filtered_asv$bin
  
  #Obtaining the filtered ASV table.
  asv0_table <- asv %>%
    select(all_of(filtered_asv$keep))
  #Calculating K0
  k_0 <- K_calculator(asv0_table)
  
  #Calculating BC and J for filtered table
  BC_0 <- mean(na.omit(bray_curtis(asv0_table)))
  J_0 <- mean(na.omit(evenness_K(asv0_table)))
  
  #Creating a list to save the values for K, BC and J
  k_values <- c(k_0)
  BC_values <- c(BC_0)
  J_values <- c(J_0)
  added_asv <- list("initial")
  #Adding a counter to save the number of the iteration
  iteration <- c(0)
  counter <- 0
  
  if (length(discarded_asv) > 0) {
    all_combinations <- unlist(
      lapply(1:length(discarded_asv), function(x) combn(discarded_asv, x, simplify = FALSE)),
      recursive = FALSE
    )
    
    #We iterate all the possible combinations
    for (combi in all_combinations) {
      counter <- counter + 1 #I keep track of iteration number
      asv_table_i <- cbind(asv0_table, asv[, combi, drop = FALSE]) #Appending the discarded asv
      
      if (any(rowSums(asv_table_i) == 0)) {
        next  # Skip the combinations if the rows sum 0, as that will f. the K calculation
      }
      
      asv_table_i <- asv_table_i[rowSums(asv_table_i) > 0, ]
      
      #Calculating K and appending to the k_values vector
      k_i <- K_calculator(asv_table_i)
      k_values <- append(k_values, k_i)
      
      ## New BC calculation
      bc_asv_i <- mean(bray_curtis(asv_table_i))
      BC_values <- append(BC_values, bc_asv_i)
      
      ## New Evennness Calculation
      j_asv_i <- mean(evenness_K(asv_table_i))
      J_values <- append(J_values, j_asv_i)
      
      #Adding the name of the combination to the added_asv list
      added_asv <- append(added_asv, list(combi))
      #Counting the iteration step
      iteration <- append(iteration, counter)
    }
    
    final_table <- asv_table_i #Saving the last table as final table.
  } else {
    final_table <- asv0_table
    print("No ASV were filtered!")
  }
  
  asv_combination_chr <- sapply(added_asv, function(x) paste(x, collapse = ","))
  
  df_results <- data.frame(
    iteration = iteration,
    asv_combination = asv_combination_chr,
    K_values = k_values,
    BC_values = BC_values,
    J_values = J_values
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    discarded_asv = discarded_asv,
    values_filtered = data.frame(K_values = k_0, BC_values = BC_0, J_values = J_0),
    values_unfiltered = data.frame(K_values = k_unfiltered, BC_values = BC_unfiltered, J_values = J_unfiltered),
    df_results = df_results,
    final_table = final_table,
    time.take = time.take
  ))
}

# Finetune Combi II: Adapted-Combinatorial
k_finetune_combi_ii <- function(ASV) {
  start.time <- Sys.time()
  
  asv <- ASV %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  # Unfiltered metrics
  k_unfiltered <- K_calculator(asv)
  BC_unfiltered <- mean(na.omit(bray_curtis(asv)))
  J_unfiltered <- mean(na.omit(evenness_K(asv)))
  
  # Filtered ASV table
  filtered_asv <- filtering_ASV(ASV)
  discarded_asv <- filtered_asv$bin
  asv0_table <- asv %>% select(all_of(filtered_asv$keep))
  
  # Metrics after filtering
  k_0 <- K_calculator(asv0_table)
  BC_0 <- mean(na.omit(bray_curtis(asv0_table)))
  J_0 <- mean(na.omit(evenness_K(asv0_table)))
  
  # Initialise tracking
  iteration <- 0
  k_values <- k_0
  BC_values <- BC_0
  J_values <- J_0
  added_asv <- "initial"
  
  best_k <- k_0
  best_combination <- "initial"
  
  # Only test combinations if something was discarded
  if (length(discarded_asv) > 0) {
    all_combinations <- unlist(
      lapply(1:min(3, length(discarded_asv)), function(x) combn(discarded_asv, x, simplify = FALSE)),
      recursive = FALSE
    )
    
    for (combi in all_combinations) {
      iteration <- iteration + 1
      asv_table_i <- cbind(asv0_table, asv[, combi, drop = FALSE])
      
      # Avoid zero-sum rows
      asv_table_i <- asv_table_i[rowSums(asv_table_i) > 0, ]
      if (nrow(asv_table_i) == 0) next
      
      k_i <- K_calculator(asv_table_i)
      bc_i <- mean(na.omit(bray_curtis(asv_table_i)))
      j_i <- mean(na.omit(evenness_K(asv_table_i)))
      
      k_values <- c(k_values, k_i)
      BC_values <- c(BC_values, bc_i)
      J_values <- c(J_values, j_i)
      added_asv <- c(added_asv, paste(combi, collapse = ","))
      
      if (!is.na(k_i) && k_i > best_k) {
        best_k <- k_i
        best_combination <- paste(combi, collapse = ",")
      }
    }
  } else {
    message("No ASVs were filtered.")
  }
  
  # Build result dataframe
  df_results <- data.frame(
    iteration = seq(0, length(k_values)-1),
    asv_combination = added_asv,
    K_values = k_values,
    BC_values = BC_values,
    J_values = J_values
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    discarded_asv = discarded_asv,
    values_filtered = data.frame(K_values = k_0, BC_values = BC_0, J_values = J_0),
    values_unfiltered = data.frame(K_values = k_unfiltered, BC_values = BC_unfiltered, J_values = J_unfiltered),
    df_results = df_results,
    best_K = best_k,
    best_combination = best_combination,
    time.take = time.take
  ))
}


# K finetune wipeout
k_finetune_wipeout<-function(ASV){
  start.time<-Sys.time()
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  
  #Creating a counter for keeping track of the iterations
  iteration<-c()
  counter<-0
  
  #Create a vector to save the name of the re-introduced ASVs in each iteration
  added_asv<-c()
  added_asv<-append(added_asv, "None_added")
  
  #Removing all columns but the firts
  asv_initial<-asv[,1,drop=FALSE]
  asv_filt_iterative <- asv_initial

  
  #Calculate the value of this "table" 
  k_0<-K_calculator(asv_initial)
  k_values<-c()
  #Creating a vector to save the K values
  k_raw<-K_calculator(asv)
  BC_raw<-mean(na.omit(bray_curtis(asv)))
  J_raw<-mean(na.omit(evenness_K(asv)))
  
  k_values<-append(k_values, k_0)
  
  BC_0<-mean(na.omit(bray_curtis(asv_initial)))
  BC_values<-c()  
  BC_values<-append(BC_values, BC_0)
  
  J_0<-mean(na.omit(evenness_K(asv_initial)))
  J_values<-c()
  J_values<-append(J_values,J_0)
  
  
  samples_id<-rownames(asv_filt_iterative)
  
  #Add this calculation as a count into the iteration
  iteration<-append(iteration, 0)
  
  for (asv_name in names(asv)[-1]){
    counter<-counter+1
    
    reint_asv<-asv[samples_id, asv_name, drop=FALSE]
    asv_table_i<-cbind(asv_filt_iterative, reint_asv)
    asv_filt_iterative<-asv_table_i
    samples_id<-rownames(asv_filt_iterative[samples_id, ,drop=FALSE])
    
    #Calculating K for the new table
    k_i<-K_calculator(asv_filt_iterative)
    BC_i<-mean(na.omit(bray_curtis(asv_filt_iterative)))
    J_i<-mean(na.omit(evenness_K(asv_filt_iterative)))
    
    added_asv<-append(added_asv, asv_name)    
    k_values<-append(k_values, k_i)
    iteration<-append(iteration, counter)
    BC_values<-append(BC_values, BC_i)
    J_values<-append(J_values, J_i)
  }
  end.time<-Sys.time()
  time.take=round(end.time-start.time,2)
  return(list(
    values_unfiltered=data.frame(K_values=k_raw, BC_values=BC_raw, J_values=J_raw),
    df_results=data.frame(iteration=iteration, added_asv=added_asv, K_values=k_values, BC_values=BC_values, J_values=J_values),
    final_table = asv_filt_iterative,
    time.take=time.take
    
  ))
}

# Add-One

k_finetune_sequential<-function(ASV){
  start.time<-Sys.time()
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  
  #Calculate unfiltered table K; BC and J
  k_unfiltered<-K_calculator(asv)
  BC_unfiltered<-mean(na.omit(bray_curtis(asv)))
  J_unfiltered<-mean(na.omit(evenness_K(asv)))
  
  #Create vectors for saving the calculated BC, K and J 
  added_asv<-c()
  K_values<-c()
  BC_values<-c()
  J_values<-c()
  
  #Filtering the table 
  filter_output<-filtering_ASV(asv)
  
  discarded_asv<-filter_output$bin
  keep_asv<-filter_output$keep
  
  #Defining filtered table
  asv0_table<-asv%>%
    select(all_of(keep_asv))
  sample_id<-rownames(asv0_table)
  
  #Calculating K, BC and J for filtered table
  k_0<-K_calculator(asv0_table)
  BC_0<-mean(na.omit(bray_curtis(asv0_table)))
  J_0<-mean(na.omit(evenness_K(asv0_table)))
  
  #Pipeline 
  for (asv_name in discarded_asv){
    asv_table_i<-cbind(asv0_table, asv[sample_id, asv_name, drop=FALSE])
    
    k_i<-K_calculator(asv_table_i)
    bc_i<-mean(na.omit(bray_curtis(asv_table_i)))
    j_i<-mean(na.omit(evenness_K(asv_table_i)))
    
    #Append 
    K_values<-append(K_values, k_i)
    BC_values<-append(BC_values, bc_i)
    J_values<-append(J_values, j_i)
    
    added_asv<-append(added_asv, asv_name)
    
  }
  end.time<-Sys.time()
  time.take<-round(end.time-start.time,2)
  return(list(
    df_results=data.frame(added_asv=added_asv, K_values=K_values , BC_values=BC_values, J_values=J_values),
    filtered_values=data.frame(K_values=k_0, BC_values=BC_0, J_values=J_0),
    unfiltered_values=data.frame(K_values=k_unfiltered, BC_values=BC_unfiltered, J_values=J_unfiltered),
    time.take=time.take
    
  ))
  
}

# Pipeline V: Semi-wipeout

filter_conservative<-function(ASV){
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  
  min_presence<-0.5
  n_samples<-nrow(asv)
  n_asv<-ncol(asv)
  
  
  keep<-c()
  bin<-c()
  
  for (i in seq_along(asv)) {
    non_zero_count <- sum(asv[[i]] != 0)
    if (non_zero_count > n_samples * min_presence) {
      keep <- c(keep, colnames(asv)[i])
    } else {
      bin <- c(bin, colnames(asv)[i])
    }
  }
  
  return(list(keep = keep, bin = bin))
}

k_finetune_semi_wipe <- function(ASV) {
  start.time <- Sys.time()
  
  asv <- ASV %>%
    mutate_if(is.integer, as.numeric) %>%
    select(where(is.numeric))
  
  #Raw table K, BC, J
  unfiltered_values <- data.frame(
    K_values = K_calculator(asv),
    BC_values = mean(na.omit(bray_curtis(asv))),
    J_values = mean(na.omit(evenness_K(asv)))
  )
  
  #Conservative filter
  asv_filter <- filter_conservative(asv)
  discarded_asv <- asv_filter$bin
  
  asv0_table <- asv %>%
    select(all_of(asv_filter$keep))

  sample_id <- rownames(asv0_table)
  
  #Filtered table K, BC, J
  filtered_values <- data.frame(
    K_values = K_calculator(asv0_table),
    BC_values = mean(na.omit(bray_curtis(asv0_table))),
    J_values = mean(na.omit(evenness_K(asv0_table)))
  )
  
  
  added_asv <- c()
  K_values <- c()
  BC_values <- c()
  J_values <- c()
  
  asv_filt_iterative <- asv0_table
  
  for (asv_name in discarded_asv) {
    asv_table_i <- cbind(asv_filt_iterative, asv[sample_id, asv_name, drop = FALSE])
    asv_filt_iterative <- asv_table_i

    K_values <- append(K_values, K_calculator(asv_table_i))
    BC_values <- append(BC_values, mean(na.omit(bray_curtis(asv_table_i))))
    J_values <- append(J_values, mean(na.omit(evenness_K(asv_table_i))))
    added_asv <- append(added_asv, asv_name)
  }
  
  iterative_df <- data.frame(
    added_asv = added_asv,
    K_values = K_values,
    BC_values = BC_values,
    J_values = J_values
  )
  
  end.time <- Sys.time()
  time.take <- round(end.time - start.time, 2)
  
  return(list(
    values_unfiltered = unfiltered_values,
    values_filtered = filtered_values,
    iterative_df = iterative_df,
    time.take = time.take
  ))
}

# k_finetune_weighted
k_finetune_weighted<-function(ASV){
  start.time<-Sys.time()
  
  #Sorting the table
  asv<-ASV%>%
    mutate_if(is.integer, as.numeric)%>%
    select(where(is.numeric))
  
  #Calculate K; BC and Eveness for unfiltered
  values_unfiltered<-data.frame(k_unfiltered=K_calculator(asv),
                                BC_unfiltered=mean(na.omit(bray_curtis(asv))),
                                J_unfiltered=mean(na.omit(evenness_K(asv))))
  
  
  #Filter the ASV table
  filter_asv<-filtering_ASV(ASV)
  
  #saving the discarded taxa id and the saved
  keep_asv<-filter_asv$keep
  discarded_asv<-filter_asv$bin
  
  #Create filtered table
  asv0_table<-asv%>%
    select(all_of(keep_asv))
  sample_id<-rownames(asv0_table)
  
  #Calculate K, BC, J for filtered table
  values_0<-data.frame(k_0=K_calculator(asv0_table),
                       bc_0=mean(na.omit(bray_curtis(asv0_table))),
                       J_0=mean(na.omit(evenness_K(asv0_table)))
  )
  #Defining the empty lists to save the fixed alpha/beta=1 results
  k_values_alpha1<-list()
  bc_values_alpha1<-list()
  j_values_alpha1<-list()
  
  k_values_beta1<-list()
  bc_values_beta1<-list()
  j_values_beta1<-list()
  
  # Calculations for fixed alpha/beta=1
  for (name_asv in discarded_asv){
    ### Fix alpha and vary beta
    alpha<-1
    beta<-seq(1,2, by=0.1)
    asv_table_i<-cbind(asv0_table, asv[sample_id, name_asv, drop=FALSE])
    asv_iterative<-asv_table_i
    
    ##Calculations BC;K and J
    k_i<-K_weighted_function(asv_iterative, alpha=alpha, beta=beta)
    k_values_alpha1[[name_asv]]<-k_i
    
    bc_i<-mean(na.omit(bray_curtis(asv_iterative)))
    bc_values_alpha1[[name_asv]]<-bc_i
    
    j_i<-mean(na.omit(evenness_K(asv_iterative)))
    j_values_alpha1[[name_asv]]<-j_i
  }
  
  for (name_asv in discarded_asv){
    ## Fix beta and vary alpha
    beta<-1
    alpha<-seq(1,2, by=0.1)
    
    asv_table_i<-cbind(asv0_table, asv[sample_id, name_asv, drop=FALSE])
    asv_iterative<-asv_table_i
    
    #calculating BC; K and J 
    k_i<-K_weighted_function(asv_iterative, alpha=alpha, beta=beta)
    k_values_beta1[[name_asv]]<-k_i
    
    bc_i<-mean(na.omit(bray_curtis(asv_iterative)))
    bc_values_beta1[[name_asv]]<-bc_i   
    
    j_i<-mean(na.omit(evenness_K(asv_iterative)))
    j_values_beta1[[name_asv]]<-j_i
    # savin the results into the list
    
  }
  
  #Calculation for alpha/beta 0.5
  k_alpha05<-list()
  k_beta05<-list()
  
  bc_alpha05<-list()
  j_alpha05<-list()
  
  bc_beta05<-list()
  j_beta05<-list()
  for (name_asv in discarded_asv){
    #Fixing alpha 
    alpha<-0.5
    beta<-seq(1,2, 0.1)
    
    asv_table_a<-cbind(asv0_table,asv[sample_id, name_asv, drop=FALSE])
    asv_table_a_it<-asv_table_a
    
    k05_a<-K_weighted_function(asv_table_a_it, alpha, beta)
    k_alpha05[[name_asv]]<-k05_a
    bc05_a<-mean(na.omit(bray_curtis(asv_table_a)))
    bc_alpha05[[name_asv]]<-bc05_a
    j05_a<-mean(na.omit(evenness_K(asv_table_a_it)))
    j_alpha05[[name_asv]]<-j05_a
    
    #Fixin beta
    alpha<-seq(1,2, by=0.1)
    beta<-0.5
    
    k05_b<-K_weighted_function(asv_table_a,alpha, beta)
    k_beta05[[name_asv]]<-k05_b    
    bc05_b<-mean(na.omit(bray_curtis(asv_table_a)))
    bc_beta05[[name_asv]]<-bc05_b    
    j05_b<-mean(na.omit(evenness_K(asv_table_a_it)))
    j_beta05[[name_asv]]<-j05_b
    
    
  }
  
  #Saving the k results per ASV for fixed alpha=1
  df_alpha1 <- as.data.frame(
    sapply(k_values_alpha1, function(x) x$k_values))
  
  #Saving the k results per ASV for fixed beta=1
  df_beta1<-as.data.frame(
    sapply(k_values_beta1, function(x) x$k_values)
  )
  
  df_alpha05<-as.data.frame(
    sapply(k_alpha05, function(x) x$k_values)
  )
  df_beta05<-as.data.frame(
    sapply(k_beta05, function(x) x$k_values)
  )
  bi<-k_values_alpha1[[1]]$b_i#Saving b values
  ai<-k_values_beta1[[1]]$a_i #Saving the a values
  
  df_alpha1<-cbind(bi, df_alpha1)
  df_beta1<-cbind(ai, df_beta1)
  
  bi<-k_alpha05[[1]]$b_i
  ai<-k_beta05[[1]]$a_i
  
  df_alpha05<-cbind(bi,df_alpha05)
  df_beta05<-cbind(ai,df_beta05)
  
  end.time=Sys.time()
  take.time=round(end.time-start.time,2)
  return(list(
    values_unfiltered=values_unfiltered,
    values_0=values_0,
    df_alpha1=df_alpha1,
    df_beta1=df_beta1,
    df_beta05=df_beta05,
    df_alpha05=df_alpha05,
    
    bc_values_alpha1=bc_values_alpha1,
    j_values_alpha1=j_values_alpha1,
    
    bc_values_beta1=bc_values_beta1,
    j_values_beta1=j_values_beta1,
    
    take.time=take.time
  ))
}


#Plots
plot_K<-function(dataframe, k_unfilt, K_filt){
  ggplot(data=dataframe, aes(added_asv, K_values))+
    geom_point(size=3, colour="#1D9AA1")+ geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=K_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=k_unfilt, linewidth=1, linetype=2, colour="#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}

plot_BC<-function(dataframe, bc_unfilt, bc_filt){
  ggplot(data=dataframe, aes(added_asv, BC_values))+
    geom_point(size=3, colour="#1D9AA1")+ geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=bc_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=bc_unfilt, linewidth=1, linetype=2, colour="#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}
plot_J <- function(dataframe, j_unfilt, j_filt){
  ggplot(data = dataframe, aes(added_asv, J_values)) +
    geom_point(size = 3, colour = "#1D9AA1") +
    geom_line(aes(group = 1), linewidth = 1.25, colour = "#1D9AA1") +
    geom_hline(yintercept = j_filt, linewidth = 1, linetype = 2, colour = "#124163") +
    geom_hline(yintercept = j_unfilt, linewidth = 1, linetype = 2, colour = "#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}


plot_K2<-function(dataframe, k_unfilt){
  ggplot(data=dataframe, aes(iteration, K_values))+
    geom_point(size=3, colour="#1D9AA1")+ geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=k_unfilt, linewidth=1, linetype=2, colour="#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}

plot_BC2<-function(dataframe, bc_unfilt){
  ggplot(data=dataframe, aes(iteration, BC_values))+
    geom_point(size=3, colour="#1D9AA1")+ geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=bc_unfilt, linewidth=1, linetype=2, colour="#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}

plot_J2<-function(dataframe, j_unfilt){
  ggplot(data=dataframe, aes(iteration, J_values))+
    geom_point(size=3, colour="#1D9AA1")+ geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=j_unfilt, linewidth=1, linetype=2, colour="#3E8853") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
}

# Cool plots for combi
cool_plot_combi_k<-function(df, k_unfilt, k_filt){
  df <- df %>%
    arrange(iteration) %>%
    mutate(asv_list = strsplit(asv_combination, ","),
           asv_combination = factor(asv_combination, levels = unique(asv_combination)))
  
  plot<-ggplot(df, aes(x = asv_list, y = K_values)) +
    geom_point(size=3, colour = "#1D9AA1") +
    geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=k_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=k_unfilt, linewidth=1, linetype=2, colour="#3E8853")+
    scale_x_upset() +
    labs(x = "ASV Combinations", y = "K value") +
    theme(base_size = 14)+
    theme_combmatrix(combmatrix.panel.point.color.fill = "#1D9AA1",
                     combmatrix.panel.line.color = "#1D9AA1",
                     combmatrix.label.make_space = FALSE)
  return (plot)
  
}



cool_plot_combi_bc<-function(df, bc_unfilt, bc_filt){
  df <- df %>%
    arrange(asv_combination) %>%
    mutate(asv_list = strsplit(asv_combination, ","),
           asv_combination = factor(asv_combination, levels = unique(asv_combination)))
  
  plot<-ggplot(df, aes(x = asv_list, y = BC_values)) +
    geom_point(size=3, colour = "#1D9AA1") +
    geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=bc_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=bc_unfilt, linewidth=1, linetype=2, colour="#3E8853")+
    scale_x_upset() +
    labs(x = "ASV Combinations", y = "Bray-Curtis value") +
    theme(base_size = 14)+
    theme_combmatrix(combmatrix.panel.point.color.fill = "#1D9AA1",
                     combmatrix.panel.line.color = "#1D9AA1",
                     combmatrix.label.make_space = FALSE)
  return (plot)
  
}

cool_plot_combi_j<-function(df, j_unfilt, j_filt){
  df <- df %>%
    arrange(asv_combination) %>%
    mutate(asv_list = strsplit(asv_combination, ","),
           asv_combination = factor(asv_combination, levels = unique(asv_combination)))
  
  plot<-ggplot(df, aes(x = asv_list, y = J_values)) +
    geom_point(size=3, colour = "#1D9AA1") +
    geom_line(aes(group=1), linewidth=1.25, colour="#1D9AA1")+
    geom_hline(yintercept=j_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=j_unfilt, linewidth=1, linetype=2, colour="#3E8853")+
    scale_x_upset() +
    labs(x = "ASV Combinations", y = "J value") +
    theme(base_size = 14)+
    theme_combmatrix(combmatrix.panel.point.color.fill = "#1D9AA1",
                     combmatrix.panel.line.color = "#1D9AA1",
                     combmatrix.label.make_space = FALSE)
  return (plot)
  
}

# Plot comparison
plot_comparison<-function(dataframe){
  ggplot(dataframe, aes(x = build_names2, y = metric_value, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c(
      "values_max" = "#1D9AA1", 
      "values_maxbc" = "#124163", 
      "values_maxj" = "#3E8853"  
    )) +
    labs(x = "Build", y = "Metric Value")
}



longer_output_k6<-function(df){
  df<-df%>%
    rename(parameter = any_of(c("ai", "bi")))
  df<-pivot_longer(df, cols=starts_with("ASV"), names_to="ASV", values_to="k_value")
  return(df)
}

longer_output_k6_b<-function(df){
  df<-df%>%
    rename(parameter = any_of(c("ai", "bi")))
  df<-pivot_longer(df, cols=starts_with("X"), names_to="ASV", values_to="k_value")
  return(df)
}

plot_vi <- function(df) {
  # Check which column exists: bi or ai
  if ("bi" %in% colnames(df)) {
    x_col <- "bi"
  } else if ("ai" %in% colnames(df)) {
    x_col <- "ai"
  } else {
    stop("Neither 'bi' nor 'ai' column found in the dataframe.")
  }
  
  if ("ASV" %in% colnames(df)) {
    z_col <- "ASV"
  } else if ("X" %in% colnames(df)) {
    z_col <- "X"
  } else {
    stop("Neither 'X' nor 'ASV' column found in the dataframe.")
  }
  # Use aes_string to map the correct column name dynamically
  plot<-ggplot(df, aes_string(x = x_col, y = "k_value", colour = z_col)) +
    geom_point(position = position_dodge(0.1), size = 2) +
    geom_path(aes(group = z_col))
  return(plot)
}

# Plots for pipeline 6 outputs

weighted_plot<-function(df){
  plot<-ggplot(df, aes(x = parameter, y = k_value, colour = condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ ASV, scales = "free_y") +
  labs(x = "Parameter value", y = "K value", colour = "Condition") +
    scale_colour_manual(values = c(
      "a1" = "#1D9AA1",
      "b1" = "#3E8853",
      "a05" = "#124163",
      "b05" = "#6EC5B8"
    )) +
  theme(base_size = 14) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  return(plot)
}
weighted_max_plot<-function(df){
  ggplot(df, aes(x = parameter, y = k_value, colour = condition)) +
    geom_line(linewidth = 1.25) +
    geom_point(size = 3) +
    facet_wrap(~ condition, scales = "free_y") +
    labs(x = "Parameter value", y = "K value", colour = "Condition") +
    scale_colour_manual(values = c(
      "a1" = "#1D9AA1",
      "b1" = "#3E8853",
      "a05" = "#124163",
      "b05" = "#6EC5B8"
    )) +
    theme(base_size = 14) +
    theme(
      strip.text = element_text(size = 10)
    )
  
}

# Plots combined 
plot_comparative_K<-function(df, K_unfilt, K_filt){
  plot<-ggplot(df, aes(x=iteration, y=K_values, fill=experiment, shape=experiment))+
    geom_point(aes(colour=experiment), size=3)+
    scale_colour_manual(values = c(
      "mbqc_or" = "#1D9AA1",
      "arranged" = "#6EC5B8"
    ))+
    geom_line(aes(group=experiment, colour=experiment), linewidth=1.25)+
    geom_hline(yintercept=K_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=K_unfilt, linewidth=1, linetype=2, colour="#3E8853")+theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
  return(plot)
}
plot_comparative_BC<-function(df, BC_unfilt, BC_filt){
  plot<-ggplot(df, aes(x=iteration, y=BC_values, fill=experiment, shape=experiment))+
    geom_point(aes(colour=experiment), size=3)+
    scale_colour_manual(values = c(
      "mbqc_or" = "#1D9AA1",
      "arranged" = "#6EC5B8"
    ))+
    geom_line(aes(group=experiment, colour=experiment), linewidth=1.25)+
    geom_hline(yintercept=BC_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=BC_unfilt, linewidth=1, linetype=2, colour="#3E8853")+theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
  return(plot)
}

plot_comparative_J<-function(df, J_unfilt, J_filt){
  plot<-ggplot(df, aes(x=iteration, y=J_values, fill=experiment, shape=experiment))+
    geom_point(aes(colour=experiment), size=3)+
    scale_colour_manual(values = c(
      "mbqc_or" = "#1D9AA1",
      "arranged" = "#6EC5B8"
    ))+
    geom_line(aes(group=experiment, colour=experiment), linewidth=1.25)+
    geom_hline(yintercept=J_filt, linewidth=1, linetype=2, colour="#124163")+
    geom_hline(yintercept=J_unfilt, linewidth=1, linetype=2, colour="#3E8853")+theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11), 
      axis.text.y = element_text(size = 12), 
      axis.title = element_text(size = 14)  
    )
  return(plot)
}


