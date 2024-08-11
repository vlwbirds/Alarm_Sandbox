####################################################################
######################## Alarm % Sandbox ###########################
####################################################################

hwk_data <- read.csv(here("hwk2.csv"))
head(hwk_data)
view(hwk_data)

##### Alarm % Pre_Height ######

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk_data  

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Calculate the unique integer "Pre_Height" values within hwk1_data
unique_pre_heights <- unique(round(hwk1_data$Pre_Height))

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_pre_heights))
sample_sizes <- numeric(length(unique_pre_heights))

# Calculate the percentage and sample size for each unique integer "Pre_Height" value within hwk1_data
for (i in 1:length(unique_pre_heights)) {
  pre_height_value <- unique_pre_heights[i]
  total_count <- sum(round(hwk1_data$Pre_Height) == pre_height_value)
  alarm_count <- sum(round(hwk1_data$Pre_Height) == pre_height_value & hwk1_data$Mid_Vocal == "alarm")
  percentage_values[i] <- (alarm_count / total_count) * 100
  sample_sizes[i] <- total_count
}

# Create a data frame for plotting
plot_data <- data.frame(Pre_Height = unique_pre_heights, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a scatter plot with Pre_Height on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Pre_Height, y = Percentage)) +
  geom_point() +  # Scatter plot
  labs(x = "Pre_Height", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) + 
  scale_x_continuous(breaks = unique_pre_heights) +  # Add integer labels to the x-axis
  
  # Add a linear regression line while keeping the line between points
  geom_smooth(method = "lm", se = FALSE, color = "blue")

p




####### Alarm% Pre_Dense #########

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk_data" CSV file
df <- hwk_data

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Calculate the unique integer "Pre_Dense" values within hwk1_data
unique_pre_dense <- unique(hwk1_data$Pre_Dense)

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_pre_dense))
sample_sizes <- numeric(length(unique_pre_dense))

# Calculate the percentage and sample size for each unique integer "Pre_Dense" value within hwk1_data
for (i in 1:length(unique_pre_dense)) {
  pre_dense_value <- unique_pre_dense[i]
  total_count <- sum(hwk1_data$Pre_Dense == pre_dense_value)
  
  # Check if total_count is zero to avoid division by zero
  if (total_count == 0) {
    percentage_values[i] <- 0  # Set percentage to 0
    sample_sizes[i] <- 0  # Set sample size to 0
  } else {
    alarm_count <- sum(hwk1_data$Pre_Dense == pre_dense_value & hwk1_data$Mid_Vocal == "alarm")
    percentage_values[i] <- (alarm_count / total_count) * 100
    sample_sizes[i] <- total_count
  }
}

# Create a data frame for plotting
plot_data <- data.frame(Pre_Dense = unique_pre_dense, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a line plot with Pre_Dense on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Pre_Dense, y = Percentage)) +
  geom_line() +
  labs(x = "Pre_Dense", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal()

# Add text labels for sample size (n) next to each point on the plot
p + geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) + 
  scale_x_continuous(breaks = unique_pre_dense)  # Add integer labels to the x-axis


### Alarm % Pre Dense Scatter Linear Regression ####
# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk_data" CSV file
df <- hwk_data

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Calculate the unique integer "Pre_Dense" values within hwk1_data
unique_pre_dense <- unique(hwk1_data$Pre_Dense)

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_pre_dense))
sample_sizes <- numeric(length(unique_pre_dense))

# Calculate the percentage and sample size for each unique integer "Pre_Dense" value within hwk1_data
for (i in 1:length(unique_pre_dense)) {
  pre_dense_value <- unique_pre_dense[i]
  total_count <- sum(hwk1_data$Pre_Dense == pre_dense_value)
  
  # Check if total_count is zero to avoid division by zero
  if (total_count == 0) {
    percentage_values[i] <- 0  # Set percentage to 0
    sample_sizes[i] <- 0  # Set sample size to 0
  } else {
    alarm_count <- sum(hwk1_data$Pre_Dense == pre_dense_value & hwk1_data$Mid_Vocal == "alarm")
    percentage_values[i] <- (alarm_count / total_count) * 100
    sample_sizes[i] <- total_count
  }
}

# Create a data frame for plotting
plot_data <- data.frame(Pre_Dense = unique_pre_dense, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a scatter plot with Pre_Dense on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Pre_Dense, y = Percentage)) +
  geom_point() +  # Scatter plot
  labs(x = "Pre_Dense", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) + 
  scale_x_continuous(breaks = unique_pre_dense)  # Add integer labels to the x-axis

# Add a linear regression line
p + geom_smooth(method = "lm", se = FALSE, color = "blue")

####### Forage Tree Height ############

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Convert "Forage_Tree_height" to numeric (if not already numeric)
hwk1_data$Forage_Tree_height <- as.numeric(hwk1_data$Forage_Tree_height)

# Filter out rows with NA and negative "Forage_Tree_height" values
hwk1_data <- hwk1_data[!is.na(hwk1_data$Forage_Tree_height) & hwk1_data$Forage_Tree_height >= 0,]

# Calculate the unique integer "Forage_Tree_height" values within hwk1_data
unique_forage_tree_heights <- unique(round(hwk1_data$Forage_Tree_height, digits = 0))  # Round to 0 decimal places

# Rest of your code (calculating percentages, creating plots, etc.) remains the same


# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_forage_tree_heights))
sample_sizes <- numeric(length(unique_forage_tree_heights))

# Calculate the percentage and sample size for each unique integer "Forage_Tree_height" value within hwk1_data
for (i in 1:length(unique_forage_tree_heights)) {
  forage_tree_height_value <- unique_forage_tree_heights[i]
  total_count <- sum(round(hwk1_data$Forage_Tree_height) == forage_tree_height_value)
  alarm_count <- sum(round(hwk1_data$Forage_Tree_height) == forage_tree_height_value & hwk1_data$Mid_Vocal == "alarm")
  percentage_values[i] <- (alarm_count / total_count) * 100
  sample_sizes[i] <- total_count
}

# Create a data frame for plotting
plot_data <- data.frame(Forage_Tree_height = unique_forage_tree_heights, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a scatter plot with Forage_Tree_height on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Forage_Tree_height, y = Percentage)) +
  geom_point() +  # Scatter plot
  labs(x = "Forage_Tree_height", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) + 
  scale_x_continuous(breaks = unique_forage_tree_heights) +  # Add integer labels to the x-axis
  
  # Add a linear regression line while keeping the line between points
  geom_smooth(method = "lm", se = FALSE, color = "blue")

p


############## Flight Proximity ################

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Convert "Flight_proximity" to numeric (if not already numeric)
hwk1_data$Flight_proximity <- as.numeric(hwk1_data$Flight_proximity)

# Filter out rows with NA and negative "Flight_proximity" values
hwk1_data <- hwk1_data[!is.na(hwk1_data$Flight_proximity) & hwk1_data$Flight_proximity >= 0,]

# Calculate the unique integer "Flight_proximity" values within hwk1_data
unique_flight_proximity <- unique(round(hwk1_data$Flight_proximity, digits = 0))  # Round to 0 decimal places

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_flight_proximity))
sample_sizes <- numeric(length(unique_flight_proximity))

# Calculate the percentage and sample size for each unique integer "Flight_proximity" value within hwk1_data
for (i in 1:length(unique_flight_proximity)) {
  flight_proximity_value <- unique_flight_proximity[i]
  total_count <- sum(round(hwk1_data$Flight_proximity) == flight_proximity_value)
  alarm_count <- sum(round(hwk1_data$Flight_proximity) == flight_proximity_value & hwk1_data$Mid_Vocal == "alarm")
  percentage_values[i] <- (alarm_count / total_count) * 100
  sample_sizes[i] <- total_count
}

# Create a data frame for plotting
plot_data <- data.frame(Flight_proximity = unique_flight_proximity, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a scatter plot with Flight_proximity on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Flight_proximity, y = Percentage)) +
  geom_point() +  # Scatter plot
  labs(x = "Flight_proximity", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) + 
  scale_x_continuous(breaks = unique_flight_proximity) +  # Add integer labels to the x-axis
  
  # Add a linear regression line while keeping the line between points
  geom_smooth(method = "lm", se = FALSE, color = "blue")

p



############ Flock Size ###############

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Preprocess the "Sex" column to make it case-insensitive and map values
hwk1_data$Sex <- tolower(hwk1_data$Sex)  # Convert to lowercase
hwk1_data$Sex <- ifelse(hwk1_data$Sex == "f", "F", ifelse(hwk1_data$Sex == "m", "M", "U"))

# Calculate the unique values within the "Sex" column
unique_sex <- unique(hwk1_data$Sex)

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_sex))
sample_sizes <- numeric(length(unique_sex))

# Calculate the percentage and sample size for each unique value within the "Sex" column
for (i in 1:length(unique_sex)) {
  sex_value <- unique_sex[i]
  total_count <- sum(hwk1_data$Sex == sex_value)
  alarm_count <- sum(hwk1_data$Sex == sex_value & hwk1_data$Mid_Vocal == "alarm")
  percentage_values[i] <- (alarm_count / total_count) * 100
  sample_sizes[i] <- total_count
}

# Create a data frame for plotting
plot_data <- data.frame(Sex = unique_sex, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a bar chart with Sex on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Sex, y = Percentage, fill = Sex)) +
  geom_bar(stat = "identity") +  # Bar chart
  labs(x = "Sex", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), position = position_stack(vjust = 0.5))  # Add text labels above bars

p


######### Veg_under_height ###########3

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Filter out rows with NA values in the "Veg_under_height" column
hwk1_data <- hwk1_data[!is.na(hwk1_data$Veg_under_height),]

# Round the values in the "Veg_under_height" column to integers
hwk1_data$Veg_under_height <- round(hwk1_data$Veg_under_height)

# Calculate the unique integer values within the "Veg_under_height" column
unique_veg_under_height <- unique(hwk1_data$Veg_under_height)

# Initialize empty vectors to store the percentages and sample sizes
percentage_values <- numeric(length(unique_veg_under_height))
sample_sizes <- numeric(length(unique_veg_under_height))

# Calculate the percentage and sample size for each unique integer value within the "Veg_under_height" column
for (i in 1:length(unique_veg_under_height)) {
  veg_under_height_value <- unique_veg_under_height[i]
  total_count <- sum(hwk1_data$Veg_under_height == veg_under_height_value)
  alarm_count <- sum(hwk1_data$Veg_under_height == veg_under_height_value & hwk1_data$Mid_Vocal == "alarm")
  percentage_values[i] <- (alarm_count / total_count) * 100
  sample_sizes[i] <- total_count
}

# Create a data frame for plotting
plot_data <- data.frame(Veg_under_height = unique_veg_under_height, Percentage = percentage_values, Sample_Size = sample_sizes)

# Create a scatter plot with Veg_under_height on the x-axis and Percentage on the y-axis
p <- ggplot(plot_data, aes(x = Veg_under_height, y = Percentage)) +
  geom_point() +  # Scatter plot
  labs(x = "Veg_under_height", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal() +
  geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5) +
  
  # Add a linear regression line using method = "lm"
  geom_smooth(method = "lm", se = FALSE, color = "blue")

p


#### 10m vert Density ####

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk" CSV file
df <- hwk_data  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Create an empty list to store the percentage values for each density meter
percentage_values_list <- list()

# Calculate the percentage of "Alarm" occurrences for all density meters combined (1 through 10)
columns_to_sum <- paste0("X5mRad_0", 1:9)
columns_to_sum <- c(columns_to_sum, "X5mRad_10")  # Add the 10th column
meter_data <- hwk1_data[rowSums(hwk1_data[, columns_to_sum]) > 0, ]
total_count <- nrow(meter_data)
alarm_count <- sum(meter_data$Mid_Vocal == "alarm")
percentage_alarm <- (alarm_count / total_count) * 100
percentage_values_list[[1]] <- percentage_alarm

# Create a data frame for plotting
plot_data <- data.frame(Density_Meter = "1-10", Percentage_Alarm = unlist(percentage_values_list))

# Create a box plot or violin plot to show the distribution of percentage values
p <- ggplot(plot_data, aes(x = Density_Meter, y = Percentage_Alarm)) +
  geom_boxplot() +  # Use geom_violin() for a violin plot
  labs(x = "Density Meters 1-10", y = "Percentage of 'Alarm' Instances (Mid_Vocal)") +
  theme_minimal()

p


######## 3d Experiment ########

# Load the rgl library
library(rgl)

# Read the data from the "hwk" CSV file
df <- hwk_data  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Create an empty list to store the data points for the 3D plot
data_points <- list()

# Calculate the percentage of "Alarm" occurrences for each combination of Pre_Height and Pre_Dense
unique_pre_heights <- unique(hwk1_data$Pre_Height)
unique_pre_denses <- unique(hwk1_data$Pre_Dense)

for (pre_height in unique_pre_heights) {
  for (pre_dense in unique_pre_denses) {
    subset_data <- hwk1_data[hwk1_data$Pre_Height == pre_height & hwk1_data$Pre_Dense == pre_dense, ]
    
    total_count <- nrow(subset_data)
    alarm_count <- sum(subset_data$Mid_Vocal == "alarm")
    percentage_alarm <- (alarm_count / total_count) * 100
    
    if (!is.na(percentage_alarm)) {
      data_points[[length(data_points) + 1]] <- c(pre_height, pre_dense, percentage_alarm)
    }
  }
}

# Convert data_points to a matrix
data_matrix <- do.call(rbind, data_points)

# Create the 3D scatter plot
plot3d(data_matrix[, 1], data_matrix[, 2], data_matrix[, 3], col = heat.colors(100), size = 2)

# Add labels to the axes
axes3d(labels = c("Pre_Height", "Pre_Dense", "Percentage Alarm"))

# Add a legend
legend3d("topright", legend = "Percentage Alarm", col = heat.colors(100), size = 0.5, pch = 16)



###### 10m Exp ###### works

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk_data" CSV file
df <- hwk_data  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Identify the "X5m" columns
x5m_columns <- grep("^X5mRad_", names(hwk1_data))

# Create a function to clean and convert a column to numeric
clean_and_convert_to_numeric <- function(column) {
  numeric_values <- as.numeric(gsub("[^0-9.]", "", column))
  numeric_values[!is.na(numeric_values)]
}

# Initialize vectors to store data
x5m_avgs <- numeric(0)
alarm_percentages <- numeric(0)

# Loop through the rows and calculate averages and percentages
for (i in 1:nrow(hwk1_data)) {
  numeric_values <- clean_and_convert_to_numeric(hwk1_data[i, x5m_columns])
  if (length(numeric_values) > 0) {
    x5m_avgs <- c(x5m_avgs, mean(numeric_values))
    total_count <- 1  # Since we are working with one row
    alarm_count <- ifelse(hwk1_data[i, "Mid_Vocal"] == "alarm", 1, 0)
    alarm_percentages <- c(alarm_percentages, (alarm_count / total_count) * 100)
  }
}

# Create a data frame with the calculated values
plot_data <- data.frame(X5m_Avg = x5m_avgs, Percentage_Alarm = alarm_percentages)

# Create a linear regression model for "Percentage_Alarm" using "X5m_Avg"
lm_model <- lm(Percentage_Alarm ~ X5m_Avg, data = plot_data)

# Create a scatter plot with the linear regression line
p <- ggplot(plot_data, aes(x = X5m_Avg, y = Percentage_Alarm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  labs(x = "Average X5m", y = "Percentage of 'Alarm' Instances (Mid_Vocal)", title = "Linear Regression for Alarm Percentage") +
  theme_minimal()

print(p)

# Print the summary of the linear regression model
summary(lm_model)

##### 10m Avg Alarm + Pre_Height ########

# Load the ggplot2 library
library(ggplot2)

# Read the data from the "hwk_data" CSV file
df <- hwk_data  # Make sure the file name is correct

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Identify the "X5m" columns
x5m_columns <- grep("^X5mRad_", names(hwk1_data))

# Create a function to clean and convert a column to numeric
clean_and_convert_to_numeric <- function(column) {
  numeric_values <- as.numeric(gsub("[^0-9.]", "", column))
  numeric_values[!is.na(numeric_values)]
}

# Initialize vectors to store data
x5m_avgs <- numeric(0)
alarm_percentages <- numeric(0)
pre_heights <- numeric(0)

# Loop through the rows and calculate averages, percentages, and pre_heights
for (i in 1:nrow(hwk1_data)) {
  numeric_values <- clean_and_convert_to_numeric(hwk1_data[i, x5m_columns])
  if (length(numeric_values) > 0) {
    x5m_avgs <- c(x5m_avgs, mean(numeric_values))
    total_count <- 1  # Since we are working with one row
    alarm_count <- ifelse(hwk1_data[i, "Mid_Vocal"] == "alarm", 1, 0)
    alarm_percentages <- c(alarm_percentages, (alarm_count / total_count) * 100)
    pre_heights <- c(pre_heights, hwk1_data[i, "Pre_Height"])
  }
}

# Create a data frame with the calculated values
plot_data <- data.frame(X5m_Avg = x5m_avgs, Percentage_Alarm = alarm_percentages, Pre_Height = pre_heights)

# Create a linear regression model for "Percentage_Alarm" using "X5m_Avg"
lm_model <- lm(Percentage_Alarm ~ X5m_Avg, data = plot_data)

# Create a scatter plot with the linear regression line
p <- ggplot(plot_data, aes(x = X5m_Avg, y = Percentage_Alarm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add linear regression line
  geom_point(data = plot_data, aes(x = X5m_Avg, y = Pre_Height), color = "red", size = 2) +  # Scatter plot for Pre_Height
  labs(x = "Average X5m", y = "Percentage of 'Alarm' Instances (Mid_Vocal)", title = "Linear Regression for Alarm Percentage") +
  theme_minimal()

print(p)

# Print the summary of the linear regression model
summary(lm_model)


###### 10m Avg #########

# Load the ggplot2 library
library(ggplot2)

# Read the data from the new CSV file "hwk_10m_avg.csv"
df <- read.csv("hwk_10m_avg.csv")

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Group the data by unique values of X5mRad_Avg and calculate the percentage of "Alarm" instances for each group
library(dplyr)
result_data <- hwk1_data %>%
  group_by(X5mRad_Avg) %>%
  summarise(Percentage_Alarm = (sum(Mid_Vocal == "alarm") / n()) * 100)

# Create a scatter plot with a regression line
p <- ggplot(result_data, aes(x = X5mRad_Avg, y = Percentage_Alarm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "X5mRad_Avg", y = "Percentage of 'Alarm' Instances (Mid_Vocal)", title = "Alarm Percentage vs. X5mRad_Avg with Regression") +
  theme_minimal()

print(p)

# Load the ggplot2 library
library(ggplot2)

# Read the data from the new CSV file "hwk_10m_avg.csv"
df <- hwk_data

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Group the data by unique values of Flight_Proximity and calculate the percentage of "Alarm" instances for each group
library(dplyr)
result_data <- hwk1_data %>%
  group_by(Flight_proximity) %>%
  summarise(Percentage_Alarm = (sum(Mid_Vocal == "alarm") / n()) * 100)

# Create a scatter plot with a parabolic curve
p <- ggplot(result_data, aes(x = Flight_proximity, y = Percentage_Alarm)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "blue") +
  labs(x = "Flight_Proximity", y = "Percentage of 'Alarm' Instances (Mid_Vocal)", title = "Alarm Percentage vs. Flight_Proximity with Parabolic Curve") +
  theme_minimal()

print(p)

### Parabolic Flight_proximity ###

# Load the ggplot2 library
library(ggplot2)

# Read the data from the new CSV file "hwk_10m_avg.csv"
df <- hwk_data

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Group the data by unique values of Flight_Proximity and calculate the percentage of "Alarm" instances for each group
library(dplyr)
result_data <- hwk1_data %>%
  group_by(Flight_proximity) %>%
  summarise(Percentage_Alarm = (sum(Mid_Vocal == "alarm") / n()) * 100,
            Sample_Size = n())  # Calculate sample size for each group

# Create a scatter plot with a cubic curve
p <- ggplot(result_data, aes(x = Flight_proximity, y = Percentage_Alarm)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "blue") +
  labs(x = "Flight_Proximity", y = "Percentage of 'Alarm' Instances (Mid_Vocal)", title = "Alarm Percentage vs. Flight_Proximity with Cubic Curve") +
  theme_minimal()

# Add sample size labels to the points
p <- p + geom_text(aes(label = Sample_Size), hjust = -0.2, vjust = 0.5)

print(p)


##### Tar_Sp Vs Mid_Vocal % ###

# Load the ggplot2 library
library(ggplot2)

# Read the data from your CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Convert Mid_Vocal values to lowercase
df$Mid_Vocal <- tolower(df$Mid_Vocal)

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Calculate the percentage of 'alarm' occurrences within each Tar_Sp category
result_data <- hwk1_data %>%
  group_by(Tar_Sp2) %>%
  summarise(Percentage_Alarm = (sum(Mid_Vocal == "alarm") / n()) * 100,
            Sample_Size = n()) %>%
  ungroup()  # Ungroup the data frame

# Create a custom color palette with more pastel colors
custom_palette <- c(
  "#FFCCCC", "#FFDDAA", "#FFEEBB", "#FFDDDD",  # Light reds and pinks
  "#FFE6CC", "#FFEEDD", "#FFEECC", "#FFFFDD",  # Light oranges and yellows
  "#AAAAAA", "#222345", "#CCFFDD", "#DDFFCC"   # Light greens
  # You can add more colors as needed
)

# Create a bar chart with a single custom fill color
p <- ggplot(result_data, aes(x = Tar_Sp2, y = Percentage_Alarm)) +
  geom_bar(stat = "identity", width = 0.5, fill = custom_palette[9]) +  # Set a specific color
  geom_text(aes(label = Sample_Size), vjust = -0.5, size = 4) +
  labs(x = "Tar_Sp", y = "Alarm Percentage", title = "Alarm Percentage by Tar_Sp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
print(p)


#####Audio Reaction by Tar_Sp########

# Load the ggplot2 library
library(ggplot2)

# Read the data from your CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1
hwk1_data <- df[df$hwk == 1, ]

# Calculate the percentage of 'alarm' occurrences within each Tar_Sp category
result_data <- hwk1_data %>%
  group_by(Tar_Sp2) %>%
  summarise(Percentage_Alarm = (sum(Audio_React == 1) / n()) * 100,
            Sample_Size = n()) %>%
  ungroup()  # Ungroup the data frame

# Create a custom color palette with more pastel colors
custom_palette <- c(
  "#FFCCCC", "#FFDDAA", "#FFEEBB", "#FFDDDD",  # Light reds and pinks
  "#FFE6CC", "#FFEEDD", "#FFEECC", "#FFFFDD",  # Light oranges and yellows
  "#AAAAAA", "#232234", "#CCFFDD", "#DDFFCC"   # Light greens
  # You can add more colors as needed
)

# Create a bar chart with a single custom fill color
p <- ggplot(result_data, aes(x = Tar_Sp2, y = Percentage_Alarm)) +
  geom_bar(stat = "identity", width = 0.5, fill = custom_palette[9]) +  # Set a specific color
  geom_text(aes(label = Sample_Size), vjust = -0.5, size = 4) +
  labs(x = "Tar_Sp", y = "Audio Reaction Percent", title = "Audio Reaction by Tar_Sp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
print(p)


####### Mid_Vocal by Tar_Sp2 ###########

# Load the required libraries
library(ggplot2)
library(dplyr)

# Read your data from the CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1 and Mid_Vocal is not NA
filtered_data <- df %>%
  filter(hwk == 1, !is.na(Mid_Vocal))

# Define the order of levels for Mid_Vocal
order_levels <- c("song","quiet", "call", "alarm")

# Reorder the Mid_Vocal factor according to the defined order
filtered_data$Mid_Vocal <- factor(filtered_data$Mid_Vocal, levels = order_levels)

# Calculate the total count within each Tar_Sp2 category
total_counts <- filtered_data %>%
  group_by(Tar_Sp2) %>%
  summarise(total = n())

# Calculate the percentage of each Mid_Vocal category within each Tar_Sp2 category
percentage_data <- filtered_data %>%
  group_by(Tar_Sp2, Mid_Vocal) %>%
  summarise(Percentage = (n() / total_counts$total) * 100) %>%
  ungroup()

# Create a bar chart with stacked bars, ensuring that each category sums to 100%
p <- ggplot(percentage_data, aes(x = Tar_Sp2, y = Percentage, fill = Mid_Vocal)) +
  geom_bar(stat = "identity", position = position_fill(vjust = 0.5), width = 0.7) +
  labs(x = "Tar_Sp2", y = "Percentage (%)", title = "Percentage by Tar_Sp2") +
  scale_fill_manual(
    values = c("alarm" = "gray36", "call" = "gray99", "quiet" = "gray60", "song" = "gray80"),
    guide = guide_legend(reverse = TRUE)  # Reverse the legend order
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  coord_flip()  # Flip the coordinates for horizontal bars

print(p)

##### Mid_Vocal by Tar_Sp2 Ordered? ########

# Load the required libraries
library(ggplot2)
library(dplyr)

# Read your data from the CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1 and Mid_Vocal is not NA
filtered_data <- df %>%
  filter(hwk == 1, !is.na(Mid_Vocal))

# Define the order of levels for Mid_Vocal
order_levels <- c("song", "quiet", "call", "alarm")

# Reorder the Mid_Vocal factor according to the defined order
filtered_data$Mid_Vocal <- factor(filtered_data$Mid_Vocal, levels = order_levels)

# Calculate the total count within each Tar_Sp2 category
total_counts <- filtered_data %>%
  group_by(Tar_Sp2) %>%
  summarise(total = n())

# Calculate the percentage of each Mid_Vocal category within each Tar_Sp2 category
percentage_data <- filtered_data %>%
  group_by(Tar_Sp2, Mid_Vocal) %>%
  summarise(Percentage = (n() / total_counts$total) * 100) %>%
  ungroup()

# Calculate the total alarm percentage for each Tar_Sp2
total_alarm_percentage <- percentage_data %>%
  filter(Mid_Vocal == "alarm") %>%
  select(Tar_Sp2, Total_Alarm_Percentage = Percentage)

# Arrange the Tar_Sp2 categories by total alarm percentage
percentage_data <- percentage_data %>%
  mutate(Tar_Sp2 = factor(Tar_Sp2, levels = total_alarm_percentage$Tar_Sp2[order(total_alarm_percentage$Total_Alarm_Percentage)]))

# Create a bar chart with stacked bars, ensuring that each category sums to 100%
p <- ggplot(percentage_data, aes(x = Tar_Sp2, y = Percentage, fill = Mid_Vocal)) +
  geom_bar(stat = "identity", position = position_fill(vjust = 0.5), width = 0.7) +
  labs(x = "Tar_Sp2", y = "Percentage (%)", title = "Percentage by Tar_Sp2") +
  scale_fill_manual(
    values = c(
      "alarm" = "gray20",
      "call" = "gray40",
      "quiet" = "gray60",
      "song" = "gray80",
      "call_pattern" = "transparent"
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentages

print(p)

##### Mid_Vocal by Tar_Sp2 with hashes? #########

# Load the required libraries
library(ggplot2)
library(dplyr)

# Read your data from the CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1 and Mid_Vocal is not NA
filtered_data <- df %>%
  filter(hwk == 1, !is.na(Mid_Vocal))

# Define the order of levels for Mid_Vocal
order_levels <- c("song","quiet", "call", "alarm")

# Reorder the Mid_Vocal factor according to the defined order
filtered_data$Mid_Vocal <- factor(filtered_data$Mid_Vocal, levels = order_levels)

# Calculate the total count within each Tar_Sp2 category
total_counts <- filtered_data %>%
  group_by(Tar_Sp2) %>%
  summarise(total = n())

# Calculate the percentage of each Mid_Vocal category within each Tar_Sp2 category
percentage_data <- filtered_data %>%
  group_by(Tar_Sp2, Mid_Vocal) %>%
  summarise(Percentage = (n() / total_counts$total) * 100) %>%
  ungroup()

# Create a bar chart with stacked bars, ensuring that each category sums to 100%
p <- ggplot(percentage_data, aes(x = Tar_Sp2, y = Percentage, fill = Mid_Vocal)) +
  geom_bar(stat = "identity", position = position_fill(vjust = 0.5), width = 0.7) +
  labs(x = "Tar_Sp2", y = "Percentage (%)", title = "Percentage by Tar_Sp2") +
  scale_fill_manual(
    values = c(
      "alarm" = "gray20", 
      "call" = "gray40", 
      "quiet" = "gray60", 
      "song" = "gray80",
      "call_pattern" = "transparent"  # Define a pattern color with a transparent fill
    ),
    guide = guide_legend(reverse = TRUE)  # Reverse the legend order
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_flip()  # Flip the coordinates for horizontal bars

print(p)

###### Mid_Vocal by Tar_Sp2, side bars ########

# Load the required libraries
library(ggplot2)
library(dplyr)

# Read your data from the CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1 and Mid_Vocal is not NA
filtered_data <- df %>%
  filter(hwk == 1, !is.na(Mid_Vocal))

# Define the order of levels for Mid_Vocal
order_levels <- c("alarm", "call", "quiet", "song")

# Reorder the Mid_Vocal factor according to the defined order
filtered_data$Mid_Vocal <- factor(filtered_data$Mid_Vocal, levels = order_levels)

# Calculate the total count within each Tar_Sp2 category
total_counts <- filtered_data %>%
  group_by(Tar_Sp2) %>%
  summarize(total = n()) %>%
  ungroup()

# Calculate the percentage of each Mid_Vocal category within each Tar_Sp2 category
percentage_data <- filtered_data %>%
  group_by(Tar_Sp2, Mid_Vocal) %>%
  summarize(Percentage = (n() / total_counts$total) * 100) %>%
  ungroup()

# Filter Tar_Sp2 to only include species with "alarm" occurrences
filtered_species <- percentage_data %>%
  filter(Mid_Vocal == "alarm") %>%
  select(Tar_Sp2)

# Create a bar chart with grouped bars and adjust the percentage axis
p <- ggplot(percentage_data, aes(x = Tar_Sp2, y = Percentage, fill = Mid_Vocal)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "Tar_Sp2", y = "Percentage (%)", title = "Percentage by Tar_Sp2") +
  scale_y_continuous(limits = c(0, 100)) +  # Set y-axis limits to 0-100%
  scale_fill_manual(
    values = c("alarm" = "gray20", "call" = "gray40", "quiet" = "gray60", "song" = "gray80"),
    guide = guide_legend(reverse = TRUE)  # Reverse the legend order
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_flip()  # Flip the coordinates for horizontal bars

print(p)


# Load the required libraries
library(ggplot2)
library(dplyr)

# Read your data from the CSV file (replace "your_data.csv" with your actual file name)
df <- hwk_data

# Filter the data for rows where hwk is equal to 1 and Mid_Vocal is not NA
filtered_data <- df %>%
  filter(hwk == 1, !is.na(Mid_Vocal))

# Define the order of levels for Mid_Vocal
order_levels <- c("alarm", "call", "quiet", "song")

# Reorder the Mid_Vocal factor according to the defined order
filtered_data$Mid_Vocal <- factor(filtered_data$Mid_Vocal, levels = order_levels)

# Calculate the total count within each Tar_Sp2 category
total_counts <- filtered_data %>%
  group_by(Tar_Sp2) %>%
  summarize(total = n()) %>%
  ungroup()

# Calculate the percentage of each Mid_Vocal category within each Tar_Sp2 category
percentage_data <- filtered_data %>%
  group_by(Tar_Sp2, Mid_Vocal) %>%
  summarize(Percentage = (n() / total_counts$total) * 100) %>%
  ungroup()

# Calculate the alarm percentage within each Tar_Sp2 category
alarm_percentage <- percentage_data %>%
  filter(Mid_Vocal == "alarm") %>%
  arrange(desc(Percentage))  # Arrange in descending order by alarm percentage

# Reorder the Tar_Sp2 factor levels based on alarm percentage
percentage_data$Tar_Sp2 <- factor(percentage_data$Tar_Sp2, levels = alarm_percentage$Tar_Sp2)

# Create a bar chart with grouped bars and adjust the percentage axis
p <- ggplot(percentage_data, aes(x = Tar_Sp2, y = Percentage, fill = Mid_Vocal)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "Tar_Sp2", y = "Percentage (%)", title = "Percentage by Tar_Sp2") +
  scale_y_continuous(limits = c(0, 100)) +  # Set y-axis limits to 0-100%
  scale_fill_manual(
    values = c("alarm" = "gray20", "call" = "gray40", "quiet" = "gray60", "song" = "gray80"),
    guide = guide_legend(reverse = TRUE)  # Reverse the legend order
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  coord_flip()  # Flip the coordinates for horizontal bars

print(p)
