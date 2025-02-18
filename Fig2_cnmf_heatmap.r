project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
library(ggplot2)
library(reshape2)
#!gsutil cp $bucket/cellranger_output_directory/bladder/harmony/cellphone/final_data_heatmap005.csv .
#df = pd.read_csv('final_data_heatmap005.csv')

system(paste0("gsutil cp -r ", bucket, "/cellranger_output_directory/bladder/harmony/cellphone/final_data_heatmap005.csv ."),intern=TRUE)
df <- read.csv('final_data_heatmap005.csv')
df

str(df$ycoord)
x_coord <- df$xcoord
y_coord <- df$ycoord
value <- df$value
foo <- data.frame(x_coord, y_coord, value)
foo
matrix_data1 <- acast(foo, y_coord ~ x_coord, value.var = "value")
matrix_data1
ggplot(melt(matrix_data1), aes(Var2, Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "X Coordinate", y = "Y Coordinate", title = "CellphoneDB result HeatMap")
q<- ggplot(df, aes(x = sample_source_target, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#1E90FF", high = "red") +
  theme_classic() #theme_void() #theme_dark() #



qq<- ggplot(df, aes(x = sample_source_target, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#1E90FF", high = "red") +
  theme_classic() #theme_void() #theme_dark() #

q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
qq + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df
q<- ggplot(df, aes(x = sample_source_target, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = target, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#79baec", high = "red") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#bdd5e7", high = "red") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "#bdd5e7", high = "#1034a6") +
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
 #scale_fill_gradient2(low = "#7382df", mid= "#ffsa10", high="red", midpoint = 0.75, na.value = "grey50")+
 #v7 #scale_fill_gradient2(low = "#fe7a15", mid="#6788f0" , high="#000395", midpoint = 0.5, na.value = "grey50")+
 scale_fill_gradient2(low = "#ff7b89", mid="#6f5f90" , high="#000395", midpoint = 1, na.value = "grey50")+

#scale_fill_gradient(low = "#f26ba6", high="#6485ee",na.value = "grey50")+
  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))





q<- ggplot(df, aes(x = sample_code, y = ycoord)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
 scale_fill_gradient2(low = "#ff9190", mid="#5e72eb" , high="#000395", midpoint = 1.3, na.value = "grey50")+
 #scale_fill_gradient(low = "#ff9190",  high="#5e72eb")+

  theme_classic() 

q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))






ggplot(df, aes(x = response, y = ligand_receptor)) +
  geom_tile(aes(fill = lr_means), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

#q + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


ggplot(df, aes(x = source, y = ligand_receptor)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()

