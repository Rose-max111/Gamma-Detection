library(rhdf5)

# 设置要创建的文件夹路径
folder_path <- "data/"

# 检查文件夹是否已存在
if (!file.exists(folder_path)) {
  # 创建文件夹
  dir.create(folder_path, recursive = TRUE)
}

# deposit_test.h5文件的数据预处理
deposit_test <- h5read("deposit_test.h5", "track_info")

for (i in 0:1999) {
  data <- deposit_test[deposit_test$event_id == i, -1]
  file_name <- paste("data/", i, ".txt", sep = "")
  
  # 在data/i.txt的第一行标出event_id为i的行数
  n <- sum(deposit_test$event_id == i)
  cat(n, "\n", file = file_name, append = FALSE)

  # 从data/i.txt的第二行开始输出event_id为i的数据
  write.table(data, file = file_name, sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
}
