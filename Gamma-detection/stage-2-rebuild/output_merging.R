library(rhdf5)
library(data.table)

deposit_test <- h5read("deposit_test.h5", "track_info")
order_test <- data.frame(event_id = as.data.frame(deposit_test)["event_id"])

# event_id为671、741、1472、1647和1984的事件跑不出结果，我们采用另一种办法来获得其order
# 即直接根据deposit_test中time的大小粗略估计其order
setDT(deposit_test)
deposit_test_sorted <- as.data.frame(deposit_test[, .SD[order(time)], by = event_id])

for (i in 0:1999) {
    sum <- sum(order_test$event_id == i)
    file_name <- paste("output/", i, "_ans.txt", sep = "")
    if (file.info(file_name)$size != 0) {
        order_test[order_test$event_id == i, "order"] <- read.table(file_name)
        order_test[order_test$event_id == i, "electron_id"] <- 0:(sum - 1)
    } else {
        order_test[order_test$event_id == i, "order"] <- 0:(sum - 1)
        order_test[order_test$event_id == i, "electron_id"] <- deposit_test_sorted[deposit_test_sorted$event_id == i, "electron_id"]
    }
}

# 对order排序，生成平台需要的格式
setDT(order_test)
order_test_sorted <- as.data.frame(order_test[, .SD[order(order)], by = event_id])

if (file.exists("order_test.h5")) {
  unlink("order_test.h5")
  cat("File", "order_test.h5", "has been deleted.\n")
}
h5write(order_test_sorted, file = "order_test.h5", "electron_order")
