# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\5-Cox")
getwd()
data <- read.csv("data with prs3.csv")

# 加载必要的R包
library(survival)    # Cox回归分析
library(caret)
library(car)         # 计算VIF
library(carData)     # 计算VIF
library(broom)

# 变量类别处理
data$ACDs[data$ACDs == 2] <- 0
data$All_cause_stroke[is.na(data$All_cause_stroke)] <- 0
data$Hemorrhagic_stroke[is.na(data$Hemorrhagic_stroke)] <- 0
data$Ischemic_stroke[is.na(data$Ischemic_stroke)] <- 0
data$TIA_[is.na(data$TIA_)] <- 0
library(dplyr)
data <- data %>%
  mutate(preference = case_when(
    preference == "wine" ~ "1",
    preference == "beer" ~ "2",
    preference == "spirits" ~ "3",
    preference == "mixed" ~ "4",
    TRUE ~ preference
  ))
data <- data %>%
  mutate(frequency = case_when(
    frequency == "Once or twice a week" ~ "1",
    frequency == "Three or four times a week" ~ "2",
    frequency == "Daily or almost daily" ~ "3",
    TRUE ~ frequency
  ))
data <- data %>%
  mutate(with_food = case_when(
    with_food == "Yes" ~ "1",
    with_food == "No" ~ "2",
    with_food == "It varies" ~ "3",
    TRUE ~ with_food
  ))
data$alcohol_risk <- as.factor(data$alcohol_risk)
data$frequency <- as.factor(data$frequency)
data$preference <- as.factor(data$preference)
data$with_food <- as.factor(data$with_food)
data$sex <- as.factor(data$sex)
data$ethnicity <- as.factor(data$ethnicity)
data$education <- as.factor(data$education)
data$smoking <- as.factor(data$smoking)
data$diet_quality <- as.factor(data$diet_quality)
data$sleep_quality <- as.factor(data$sleep_quality)
data$blood_pressure <- as.factor(data$blood_pressure)
data$BMI <- as.factor(data$BMI)
data$health_rating <- as.factor(data$health_rating)
str(data)

#### 排除随访时间＜2年、健康自评poor(编码4)的人群后进行敏感性分析
# 统计随访时间＜2年的人数
less2 <- sum(data$follow_up_time < 2)
cat("follow_up_time 小于 2 的行数：", less2, "\n")
data <- data[data$follow_up_time >= 2, ]
# 统计health_rating列中值为4的人数
count_4 <- sum(data$health_rating == 4)
cat("health_rating列中值为4的人数为：", count_4, "\n")
data <- data[data$health_rating != 4, ]

##### 模型1敏感性分析

# 由于按照原始单元进行分析后得出的HR值非常小，接近于1.000
# 因此对alcohol_pw进行处理，进行缩放，使得HR对应的自变量增量为1SD
data$alcohol_pw_std <- scale(data$alcohol_pw)  

# 检查模型1协变量之间的多重共线性
pw_model_1 <- vif(lm(alcohol_pw_std ~ preference + with_food, data = data))
print(pw_model_1)

## 饮酒量作为连续型变量（alcohol_pw）
# alcohol_pw与ACDs
pw_ACDs_1 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_ACDs_1)
# 整理模型结果
library(broom)
pw_ACDs_1r <- tidy(pw_ACDs_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_ACDs_1r, "pw_ACDs_model 1.csv", row.names = FALSE)

## alcohol_pw与All_cause_stroke
pw_AS_1 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_AS_1)
pw_AS_1r <- tidy(pw_AS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_1r, "pw_AS_model 1.csv", row.names = FALSE)

## alcohol_pw与Hemorrhagic_stroke
pw_HS_1 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_HS_1)
pw_HS_1r <- tidy(pw_HS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_1r, "pw_HS_model 1.csv", row.names = FALSE)

## alcohol_pw与Ischemic_stroke
pw_IS_1 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_IS_1)
pw_IS_1r <- tidy(pw_IS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_1r, "pw_IS_model 1.csv", row.names = FALSE)

## alcohol_pw与TIA
pw_TIA_1 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_TIA_1)
pw_TIA_1r <- tidy(pw_TIA_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_1r, "pw_TIA_model 1.csv", row.names = FALSE)

## 饮酒量作为分类变量（alcohol_risk）
install.packages("caret")
library(caret)
## alcohol_risk的参照组设置为1-lower risk
data$alcohol_risk <- relevel(data$alcohol_risk, ref = "1")

# 创建哑变量（Dummy Variables），去掉截距项
dummy_model <- dummyVars(~ preference + with_food, data = data, fullRank = TRUE)
df_vif <- predict(dummy_model, newdata = data)
df_vif <- as.data.frame(df_vif)  # 转换为数据框

# 转换为数据框并确保行数一致
df_vif <- as.data.frame(df_vif)
if (nrow(df_vif) != nrow(data)) {
  stop("哑变量行数与原始数据不匹配！请检查缺失值或数据处理步骤。")
}
if (any(is.na(df_vif))) {
  stop("VIF 计算的数据框仍然存在 NA，请检查数据处理流程！")
}

# 修正变量名（替换 . 避免影响公式解析）
colnames(df_vif) <- gsub("\\.", "_", colnames(df_vif))

# 计算 VIF（线性回归不需要因变量）
vif_values <- vif(lm(I(rep(1, nrow(df_vif))) ~ ., data = df_vif))  # 使用 I(rep(1, nrow(df_vif))) 避免格式问题
print(vif_values)

# alcohol_risk与ACDs
risk_ACDs_1 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_ACDs_1)
risk_ACDs_1r <- tidy(risk_ACDs_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_ACDs_1r, "risk_ACDs_model 1.csv", row.names = FALSE)

# alcohol_risk与All_cause_stroke
risk_AS_1 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_AS_1)
risk_AS_1r <- tidy(risk_AS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_1r, "risk_AS_model 1.csv", row.names = FALSE)

# alcohol_risk与Hemorrhagic_stroke
risk_HS_1 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_HS_1)
risk_HS_1r <- tidy(risk_HS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_1r, "risk_HS_model 1.csv", row.names = FALSE)

# alcohol_risk与Ischemic_stroke
risk_IS_1 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_IS_1)
risk_IS_1r <- tidy(risk_IS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_1r, "risk_IS_model 1.csv", row.names = FALSE)

# alcohol_risk与TIA
risk_TIA_1 <- coxph(Surv(follow_up_time, TIA_) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_TIA_1)
risk_TIA_1r <- tidy(risk_TIA_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_1r, "risk_TIA_model 1.csv", row.names = FALSE)

## Cox模型2调整的协变量：饮酒模式+社会人口学因素
## 饮酒作为连续型变量（alcohol_pw）
# 检查模型2协变量之间的多重共线性
pw_model_2 <- vif(lm(alcohol_pw_std ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data))
print(pw_model_2)

# alcohol_pw与ACDs
pw_ACDs_2 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_ACDs_2)
# 整理模型结果
pw_ACDs_2r <- tidy(pw_ACDs_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_ACDs_2r, "pw_ACDs_model 2.csv", row.names = FALSE)

# alcohol_pw与All_cause_stroke
pw_AS_2 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_AS_2)
# 整理模型结果
pw_AS_2r <- tidy(pw_AS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_2r, "pw_AS_model 2.csv", row.names = FALSE)

# alcohol_pw与Hemorrhagic_stroke
pw_HS_2 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_HS_2)
# 整理模型结果
pw_HS_2r <- tidy(pw_HS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_2r, "pw_HS_model 2.csv", row.names = FALSE)

# alcohol_pw与Ischemic_stroke
pw_IS_2 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_IS_2)
# 整理模型结果
pw_IS_2r <- tidy(pw_IS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_2r, "pw_IS_model 2.csv", row.names = FALSE)

# alcohol_pw与TIA
pw_TIA_2 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_TIA_2)
# 整理模型结果
pw_TIA_2r <- tidy(pw_TIA_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_2r, "pw_TIA_model 2.csv", row.names = FALSE)

## 饮酒作为分类变量（alcohol_risk）

# 创建哑变量（Dummy Variables），去掉截距项
dummy_model_2 <- dummyVars(~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data, fullRank = TRUE)
df_vif_2 <- predict(dummy_model_2, newdata = data)
df_vif_2 <- as.data.frame(df_vif_2)  # 转换为数据框
# 转换为数据框并确保行数一致
df_vif_2 <- as.data.frame(df_vif_2)
if (nrow(df_vif_2) != nrow(data)) {
  stop("哑变量行数与原始数据不匹配！请检查缺失值或数据处理步骤。")
}
if (any(is.na(df_vif_2))) {
  stop("VIF 计算的数据框仍然存在 NA，请检查数据处理流程！")
}
# 修正变量名（替换 . 避免影响公式解析）
colnames(df_vif_2) <- gsub("\\.", "_", colnames(df_vif_2))
# 计算 VIF（线性回归不需要因变量）
vif_values_2 <- vif(lm(I(rep(1, nrow(df_vif_2))) ~ ., data = df_vif_2))  # 使用 I(rep(1, nrow(df_vif))) 避免格式问题
print(vif_values_2)

## 由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并
table(data$sleep_quality)
chisq.test(table(data$sleep_quality))
data <- data %>%
  mutate(sleep_quality_mer = case_when(
    sleep_quality %in% c(1, 2) ~ 1,
    sleep_quality == 3 ~ 2,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
data$sleep_quality_mer <- as.factor(data$sleep_quality_mer)
# 创建哑变量（Dummy Variables），去掉截距项
dummy_model_2m <- dummyVars(~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data = data, fullRank = TRUE)
df_vif_2m <- predict(dummy_model_2m, newdata = data)
df_vif_2m <- as.data.frame(df_vif_2m)  # 转换为数据框
# 转换为数据框并确保行数一致
df_vif_2m <- as.data.frame(df_vif_2m)
if (nrow(df_vif_2m) != nrow(data)) {
  stop("哑变量行数与原始数据不匹配！请检查缺失值或数据处理步骤。")
}
if (any(is.na(df_vif_2m))) {
  stop("VIF 计算的数据框仍然存在 NA，请检查数据处理流程！")
}
# 修正变量名（替换 . 避免影响公式解析）
colnames(df_vif_2m) <- gsub("\\.", "_", colnames(df_vif_2m))
# 计算 VIF（线性回归不需要因变量）
vif_values_2m <- vif(lm(I(rep(1, nrow(df_vif_2m))) ~ ., data = df_vif_2m))  # 使用 I(rep(1, nrow(df_vif))) 避免格式问题
print(vif_values_2m)

# alcohol_risk与ACDs
risk_ACDs_2 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_ACDs_2)
risk_ACDs_2r <- tidy(risk_ACDs_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_ACDs_2r, "risk_ACDs_model 2.csv", row.names = FALSE)

# alcohol_risk与All_cause_stroke
risk_AS_2 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_AS_2)
risk_AS_2r <- tidy(risk_AS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_2r, "risk_AS_model 2.csv", row.names = FALSE)

# alcohol_risk与Hemorrhagic_stroke
risk_HS_2 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_HS_2)
risk_HS_2r <- tidy(risk_HS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_2r, "risk_HS_model 2.csv", row.names = FALSE)

# alcohol_risk与Ischemic_stroke
risk_IS_2 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_IS_2)
risk_IS_2r <- tidy(risk_IS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_2r, "risk_IS_model 2.csv", row.names = FALSE)

# alcohol_risk与TIA
risk_TIA_2 <- coxph(Surv(follow_up_time, TIA_) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_TIA_2)
risk_TIA_2r <- tidy(risk_TIA_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_2r, "risk_TIA_model 2.csv", row.names = FALSE)

## Cox模型3调整的协变量：饮酒模式+社会人口学因素+疾病危险因素
# 检查模型3协变量之间的多重共线性
pw_model_3 <- vif(lm(alcohol_pw_std ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + cholesterol + triglycerides + HDL + LDL + glucose + urate , data = data))
print(pw_model_3)
# 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
pw_model_4 <- vif(lm(alcohol_pw_std ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data))
print(pw_model_4)

# alcohol_pw与ACDs
pw_ACDs_3 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_ACDs_3)
# 整理模型结果
pw_ACDs_3r <- tidy(pw_ACDs_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_ACDs_3r, "pw_ACDs_model 3.csv", row.names = FALSE)

# alcohol_pw与All_cause_stroke
pw_AS_3 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_AS_3)
# 整理模型结果
pw_AS_3r <- tidy(pw_AS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_3r, "pw_AS_model 3.csv", row.names = FALSE)

# alcohol_pw与Hemorrhagic_stroke
pw_HS_3 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_HS_3)
# 整理模型结果
pw_HS_3r <- tidy(pw_HS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_3r, "pw_HS_model 3.csv", row.names = FALSE)

# alcohol_pw与Ischemic_stroke
pw_IS_3 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_IS_3)
# 整理模型结果
pw_IS_3r <- tidy(pw_IS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_3r, "pw_IS_model 3.csv", row.names = FALSE)

# alcohol_pw与TIA
pw_TIA_3 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_TIA_3)
# 整理模型结果
pw_TIA_3r <- tidy(pw_TIA_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_3r, "pw_TIA_model 3.csv", row.names = FALSE)

## 饮酒作为分类变量（alcohol_risk）
# 创建哑变量（Dummy Variables），去掉截距项
dummy_model_3 <- dummyVars(~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + cholesterol + triglycerides + HDL + LDL + glucose + urate , data = data, fullRank = TRUE)
df_vif_3 <- predict(dummy_model_3, newdata = data)
df_vif_3 <- as.data.frame(df_vif_3)  # 转换为数据框
# 转换为数据框并确保行数一致
df_vif_3 <- as.data.frame(df_vif_3)
if (nrow(df_vif_3) != nrow(data)) {
  stop("哑变量行数与原始数据不匹配！请检查缺失值或数据处理步骤。")
}
if (any(is.na(df_vif_3))) {
  stop("VIF 计算的数据框仍然存在 NA，请检查数据处理流程！")
}
# 修正变量名（替换 . 避免影响公式解析）
colnames(df_vif_3) <- gsub("\\.", "_", colnames(df_vif_3))
# 检查每个变量的唯一值数量
sapply(df_vif_3, function(x) length(unique(x)))
# 检查是否有变量的值完全相同
apply(df_vif_3, 2, function(x) var(x) == 0)  # 如果为TRUE，表示该变量的方差为0
# 移除方差为0的变量
variance_zero <- apply(df_vif_3, 2, function(x) var(x) == 0)  # 检查哪些变量的方差为0
df_vif_3 <- df_vif_3[, !variance_zero]  # 移除方差为0的变量
# 查看移除后的变量
cat("移除方差为0的变量后，剩余变量：", names(df_vif_3), "\n")
# 重新计算VIF
model <- lm(I(rep(1, nrow(df_vif_3))) ~ ., data = df_vif_3)
vif_values_3 <- vif(model)
print(vif_values_3)

## 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
## 由于由于sleep_quality_2、sleep_quality_3的VIF大于2 需要进行成分合并（将1和2合并为1-poor or moderate，3改为2-healthy）
## blood_pressure3、blood_pressure4的VIF大于2 需要进行成分合并（将1和2合并为1-normal or elevated，3改为2-stage 1, 4改为3-stage 2）
table(data$blood_pressure)
chisq.test(table(data$blood_pressure))
data <- data %>%
  mutate(blood_pressure_mer = case_when(
    blood_pressure %in% c(1, 2) ~ 1,
    blood_pressure == 3 ~ 2,
    blood_pressure == 4 ~ 3,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
data$blood_pressure_mer <- as.factor(data$blood_pressure_mer)
chisq.test(table(data$blood_pressure_mer))


# 创建哑变量（Dummy Variables），去掉截距项
dummy_model_4 <- dummyVars(~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data, fullRank = TRUE)
df_vif_4 <- predict(dummy_model_4, newdata = data)
df_vif_4 <- as.data.frame(df_vif_4)  # 转换为数据框

# 转换为数据框并确保行数一致
df_vif_4 <- as.data.frame(df_vif_4)
if (nrow(df_vif_4) != nrow(data)) {
  stop("哑变量行数与原始数据不匹配！请检查缺失值或数据处理步骤。")
}
if (any(is.na(df_vif_4))) {
  stop("VIF 计算的数据框仍然存在 NA，请检查数据处理流程！")
}

# 修正变量名（替换 . 避免影响公式解析）
colnames(df_vif_4) <- gsub("\\.", "_", colnames(df_vif_4))
# 检查每个变量的唯一值数量
sapply(df_vif_4, function(x) length(unique(x)))
# 检查是否有变量的值完全相同
apply(df_vif_4, 2, function(x) var(x) == 0)  # 如果为TRUE，表示该变量的方差为0
# 移除方差为0的变量
variance_zero <- apply(df_vif_4, 2, function(x) var(x) == 0)  # 检查哪些变量的方差为0
df_vif_4 <- df_vif_4[, !variance_zero]  # 移除方差为0的变量
# 查看移除后的变量
cat("移除方差为0的变量后，剩余变量：", names(df_vif_4), "\n")
# 重新计算VIF
model2 <- lm(I(rep(1, nrow(df_vif_4))) ~ ., data = df_vif_4)
vif_values_4 <- vif(model2)
print(vif_values_4)

# alcohol_risk与ACDs
risk_ACDs_3 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_ACDs_3)
risk_ACDs_3r <- tidy(risk_ACDs_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_ACDs_3r, "risk_ACDs_model 3.csv", row.names = FALSE)

# alcohol_risk与All_cause_stroke
risk_AS_3 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_AS_3)
risk_AS_3r <- tidy(risk_AS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_3r, "risk_AS_model 3.csv", row.names = FALSE)

# alcohol_risk与Hemorrhagic_stroke
risk_HS_3 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_HS_3)
risk_HS_3r <- tidy(risk_HS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_3r, "risk_HS_model 3.csv", row.names = FALSE)

# alcohol_risk与Ischemic_stroke
risk_IS_3 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_IS_3)
risk_IS_3r <- tidy(risk_IS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_3r, "risk_IS_model 3.csv", row.names = FALSE)

# alcohol_risk与TIA
risk_TIA_3 <- coxph(Surv(follow_up_time, TIA_) ~ alcohol_risk + preference + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_TIA_3)
risk_TIA_3r <- tidy(risk_TIA_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_3r, "risk_TIA_model 3.csv", row.names = FALSE)
