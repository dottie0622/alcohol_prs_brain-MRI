# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\5-Cox")
getwd()

data <- read.csv("data with prs3_r.csv")


##### 饮酒与疾病之间关系的初步分析
## 疾病分为ACDs、AS、HS、IS、TIA

## Cox模型1调整的协变量：酒精偏好、是否与食物共同服用
## Cox分析之前检查协变量之间的共线性关系

# 加载必要的R包
library(survival)    # Cox回归分析
library(car)         # 计算VIF
library(carData)     # 计算VIF


# 读取数据（数据框为data）
# 变量：
#   follow_up_time: 结局事件发生的时间(生存时间),连续型变量
#   ACDs、All_cause_stroke、Hemorrhagic_stroke、Ischemic_stroke、TIA_: 结局变量(是否患病, 1=患病, 0=未患病),连续型变量
#   alcohol_risk: 饮酒量分度（自变量），分类变量
#   alcohol_pw: 饮酒量，连续型变量
#   preference: 饮酒偏好(1=wine, 2=beer, 3=spirits,4=mixed)，分类变量
#   with_food: 是否与食物同时饮酒(1=Yes,2=No,3=It varies)，分类变量

# 变量类别处理
data$ACDs[data$ACDs == 2] <- 0
data$All_cause_stroke[is.na(data$All_cause_stroke)] <- 0
data$Hemorrhagic_stroke[is.na(data$Hemorrhagic_stroke)] <- 0
data$Ischemic_stroke[is.na(data$Ischemic_stroke)] <- 0
data$TIA[is.na(data$TIA)] <- 0

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

##### 饮酒量与疾病的Cox分析：模型1

# 由于按照原始单元进行分析后得出的HR值非常小，接近于1.000
# 因此对alcohol_pw进行处理，进行缩放，使得HR对应的自变量增量为1SD
data$alcohol_pw_std <- scale(data$alcohol_pw)  


# 检查模型1协变量之间的多重共线性
pw_model_1 <- vif(lm(alcohol_pw_std ~ preference + with_food, data = data))
print(pw_model_1)

## 饮酒量作为连续型变量（alcohol_pw_std）
# alcohol_pw与ACDs
pw_ACDs_1 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_ACDs_1)
# 整理模型结果
library(broom)
pw_ACDs_1r <- tidy(pw_ACDs_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_ACDs_1r, "pw_ACDs_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_ACDs_1)
print(sch_residual_test)  # 打印检验结果

## alcohol_pw与All_cause_stroke
pw_AS_1 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_AS_1)
pw_AS_1r <- tidy(pw_AS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_1r, "pw_AS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_AS_1)
print(sch_residual_test)  # 打印检验结果

## alcohol_pw与Hemorrhagic_stroke
pw_HS_1 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_HS_1)
pw_HS_1r <- tidy(pw_HS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_1r, "pw_HS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_HS_1)
print(sch_residual_test)  # 打印检验结果

## alcohol_pw与Ischemic_stroke
pw_IS_1 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_IS_1)
pw_IS_1r <- tidy(pw_IS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_1r, "pw_IS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_IS_1)
print(sch_residual_test)  # 打印检验结果

## alcohol_pw与TIA
pw_TIA_1 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food, data = data)
summary(pw_TIA_1)
pw_TIA_1r <- tidy(pw_TIA_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_1r, "pw_TIA_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_TIA_1)
print(sch_residual_test)  # 打印检验结果


## 饮酒量作为分类变量（alcohol_risk）
install.packages("caret")
library(caret)
## alcohol_risk的参照组设置为1-low risk
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
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_ACDs_1)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与All_cause_stroke
risk_AS_1 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_AS_1)
risk_AS_1r <- tidy(risk_AS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_1r, "risk_AS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_AS_1)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Hemorrhagic_stroke
risk_HS_1 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_HS_1)
risk_HS_1r <- tidy(risk_HS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_1r, "risk_HS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_HS_1)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Ischemic_stroke
risk_IS_1 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_IS_1)
risk_IS_1r <- tidy(risk_IS_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_1r, "risk_IS_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_IS_1)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与TIA
risk_TIA_1 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food, data=data)
summary(risk_TIA_1)
risk_TIA_1r <- tidy(risk_TIA_1, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_1r, "risk_TIA_model 1.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_TIA_1)
print(sch_residual_test)  # 打印检验结果

## Cox模型2调整的协变量：饮酒模式+社会人口学因素
## 饮酒模式（preference、with_food）
## 社会人口学因素（age、sex、ethnicity、education、TDI、smoking、diet_quality、sleep_quality）

## 饮酒作为连续型变量（alcohol_pw_std）
# 检查模型2协变量之间的多重共线性
pw_model_2 <- vif(lm(alcohol_pw_std ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data))
print(pw_model_2)

# alcohol_pw与ACDs
pw_ACDs_2 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_ACDs_2)
# 整理模型结果
pw_ACDs_2r <- tidy(pw_ACDs_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_ACDs_2r, "pw_ACDs_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_ACDs_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与All_cause_stroke
pw_AS_2 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_AS_2)
# 整理模型结果
pw_AS_2r <- tidy(pw_AS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_2r, "pw_AS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_AS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与Hemorrhagic_stroke
pw_HS_2 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_HS_2)
# 整理模型结果
pw_HS_2r <- tidy(pw_HS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_2r, "pw_HS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_HS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与Ischemic_stroke
pw_IS_2 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_IS_2)
# 整理模型结果
pw_IS_2r <- tidy(pw_IS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_2r, "pw_IS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_IS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与TIA
pw_TIA_2 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = data)
summary(pw_TIA_2)
# 整理模型结果
pw_TIA_2r <- tidy(pw_TIA_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_2r, "pw_TIA_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_TIA_2)
print(sch_residual_test)  # 打印检验结果

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
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_ACDs_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与All_cause_stroke
risk_AS_2 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_AS_2)
risk_AS_2r <- tidy(risk_AS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_2r, "risk_AS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_AS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Hemorrhagic_stroke
risk_HS_2 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_HS_2)
risk_HS_2r <- tidy(risk_HS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_2r, "risk_HS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_HS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Ischemic_stroke
risk_IS_2 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_IS_2)
risk_IS_2r <- tidy(risk_IS_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_2r, "risk_IS_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_IS_2)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与TIA
risk_TIA_2 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer , data=data)
summary(risk_TIA_2)
risk_TIA_2r <- tidy(risk_TIA_2, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_2r, "risk_TIA_model 2.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_TIA_2)
print(sch_residual_test)  # 打印检验结果

## Cox模型3调整的协变量：饮酒模式+社会人口学因素+疾病危险因素
## 饮酒模式（preference、with_food）
## 社会人口学因素（age、sex、ethnicity、education、TDI、smoking、diet_quality、sleep_quality）
## 健康相关危险因素（blood_pressure、BMI、health_rating、cholesterol、triglycerides、HDL、LDL、glucose、urate）

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
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_ACDs_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与All_cause_stroke
pw_AS_3 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_AS_3)
# 整理模型结果
pw_AS_3r <- tidy(pw_AS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_AS_3r, "pw_AS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_AS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与Hemorrhagic_stroke
pw_HS_3 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_HS_3)
# 整理模型结果
pw_HS_3r <- tidy(pw_HS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_HS_3r, "pw_HS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_HS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与Ischemic_stroke
pw_IS_3 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_IS_3)
# 整理模型结果
pw_IS_3r <- tidy(pw_IS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_IS_3r, "pw_IS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_IS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_pw与TIA
pw_TIA_3 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data)
summary(pw_TIA_3)
# 整理模型结果
pw_TIA_3r <- tidy(pw_TIA_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(pw_TIA_3r, "pw_TIA_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(pw_TIA_3)
print(sch_residual_test)  # 打印检验结果

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

# 计算 VIF（线性回归不需要因变量）
vif_values_3 <- vif(lm(I(rep(1, nrow(df_vif_3))) ~ ., data = df_vif_3))  # 使用 I(rep(1, nrow(df_vif))) 避免格式问题
print(vif_values_3)

## 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
## 由于由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并（将1和2合并为1-poor or moderate，3改为2-healthy）
## blood_pressure3、blood_pressure4的VIF大于5 需要进行成分合并（将1和2合并为1-normal or elevated，3改为2-stage 1, 4改为3-stage 2）
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

# 计算 VIF（线性回归不需要因变量）
vif_values_4 <- vif(lm(I(rep(1, nrow(df_vif_4))) ~ ., data = df_vif_4))  # 使用 I(rep(1, nrow(df_vif))) 避免格式问题
print(vif_values_4)

# alcohol_risk与ACDs
risk_ACDs_3 <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_ACDs_3)
risk_ACDs_3r <- tidy(risk_ACDs_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_ACDs_3r, "risk_ACDs_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_ACDs_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与All_cause_stroke
risk_AS_3 <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_AS_3)
risk_AS_3r <- tidy(risk_AS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_AS_3r, "risk_AS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_AS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Hemorrhagic_stroke
risk_HS_3 <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_HS_3)
risk_HS_3r <- tidy(risk_HS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_HS_3r, "risk_HS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_HS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与Ischemic_stroke
risk_IS_3 <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_IS_3)
risk_IS_3r <- tidy(risk_IS_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_IS_3r, "risk_IS_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_IS_3)
print(sch_residual_test)  # 打印检验结果

# alcohol_risk与TIA
risk_TIA_3 <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data=data)
summary(risk_TIA_3)
risk_TIA_3r <- tidy(risk_TIA_3, exponentiate = TRUE, conf.int = TRUE)
write.csv(risk_TIA_3r, "risk_TIA_model 3.csv", row.names = FALSE)
# 进行Schoenfeld残差检验
sch_residual_test <- cox.zph(risk_TIA_3)
print(sch_residual_test)  # 打印检验结果

##### 非线性关系分析
##### 利用限制三次样条模型探讨饮酒量与每种疾病风险的非线性关系

# 加载所需的R包
install.packages("rms")
library(rms)    # 用于限制三次样条和Cox回归分析
library(survival)  # 用于生存分析
library(ggplot2)   # 用于绘图

# 设置数据环境
ddist <- datadist(data)  # 将数据分布信息存储到ddist中
options(datadist = "ddist")  # 设置全局选项，使模型能够使用数据分布信息
summary(data$alcohol_pw_std)
# 检查 alcohol_pw_std 的数据结构
str(data$alcohol_pw_std)
# 如果 alcohol_pw_std 是矩阵，提取第一列作为向量
if (is.matrix(data$alcohol_pw_std)) {
  data$alcohol_pw_std <- data$alcohol_pw_std[, 1]
}

# 模型3
# 构建限制三次样条模型
# ACDs: 使用 rcs() 函数对饮酒量进行限制三次样条变换，节点数设为3（可根据需要调整）
ACDs_model <- cph(
  Surv(follow_up_time, ACDs) ~ rcs(alcohol_pw, 3) + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate ,
  data = data
)
# 查看模型摘要
summary(ACDs_model)  # 输出模型的详细结果，包括各变量的系数和显著性
# 检查非线性项是否显著
anova(ACDs_model)  # 对模型进行方差分析，检查饮酒量的非线性项是否显著

# 绘制限制三次样条图
# 使用 Predict() 函数生成饮酒量与疾病风险的预测值
pred_ACDs <- Predict(ACDs_model, alcohol_pw)  # 生成预测值
# 将 Log HR 转换为 HR
pred_ACDs$yhat <- exp(pred_ACDs$yhat)  # 对 yhat 进行指数变换
pred_ACDs$lower <- exp(pred_ACDs$lower)  # 对置信区间下限进行指数变换
pred_ACDs$upper <- exp(pred_ACDs$upper)  # 对置信区间上限进行指数变换
# 找到最小风险剂量（饮酒量对应的HR最低点）
min_risk_index <- which.min(pred_ACDs$yhat)
min_risk_alcohol <- pred_ACDs$alcohol_pw[min_risk_index]
min_risk_hr <- pred_ACDs$yhat[min_risk_index]
# 找到风险显著增加的阈值（HR置信区间下限超过1的最小饮酒量）
threshold_index <- which(pred_ACDs$lower > 1)
if (length(threshold_index) > 0) {
  threshold_alcohol <- pred_ACDs$alcohol_pw[threshold_index[1]]
  threshold_hr <- pred_ACDs$yhat[threshold_index[1]]
} else {
  threshold_alcohol <- NA
  threshold_hr <- NA
}

# 输出结果
cat("最小风险剂量（饮酒量）:", min_risk_alcohol, "\n")
cat("对应的危险比（HR）:", min_risk_hr, "\n")
if (!is.na(threshold_alcohol)) {
  cat("风险显著增加的阈值（饮酒量）:", threshold_alcohol, "\n")
  cat("对应的危险比（HR）:", threshold_hr, "\n")
} else {
  cat("未找到风险显著增加的阈值（置信区间未超过1）。\n")
}
# 绘制限制三次样条图
ACDs_nonlinear <- ggplot(pred_ACDs, aes(x = alcohol_pw_std, y = yhat)) +
  geom_line() +  # 绘制曲线
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  # 添加置信区间
  geom_vline(xintercept = min_risk_alcohol, linetype = "dashed", color = "blue") +  # 最小风险剂量
  geom_vline(xintercept = threshold_alcohol, linetype = "dashed", color = "red") +  # 风险显著增加的阈值
  labs(
    title = "The relationship between alcohol intake and incident ACDs risk.",
    x = "Alcohol intake",
    y = "HR"
  ) +
  theme_minimal()  # 使用简洁的主题
# 显示图形
print(ACDs_nonlinear)
ggsave("non-linear ACDs.pdf", plot = ACDs_nonlinear, width = 8, height = 6, dpi = 300)

# AS: 使用 rcs() 函数对饮酒量进行限制三次样条变换，节点数设为3（可根据需要调整）
AS_model <- cph(
  Surv(follow_up_time, All_cause_stroke) ~ rcs(alcohol_pw_std, 3) + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate ,
  data = data
)
# 查看模型摘要
summary(AS_model)  # 输出模型的详细结果，包括各变量的系数和显著性
# 检查非线性项是否显著
anova(AS_model)  # 对模型进行方差分析，检查饮酒量的非线性项是否显著

# 绘制限制三次样条图
# 使用 Predict() 函数生成饮酒量与疾病风险的预测值
pred_AS <- Predict(AS_model, alcohol_pw_std)  # 生成预测值
# 将 Log HR 转换为 HR
pred_AS$yhat <- exp(pred_AS$yhat)  # 对 yhat 进行指数变换
pred_AS$lower <- exp(pred_AS$lower)  # 对置信区间下限进行指数变换
pred_AS$upper <- exp(pred_AS$upper)  # 对置信区间上限进行指数变换
# 保存图形到本地文件
AS_nonlinear <- ggplot(pred_AS) +  # 使用ggplot2绘制图形
  geom_line() +  # 绘制曲线
  labs(
    title = "The relationship between alcohol intake and incident All-cause stroke risk.",
    x = "Alcohol intake",
    y = "HR"
  ) +
  theme_minimal()  # 使用简洁的主题
ggsave("non-linear AS.pdf", plot = AS_nonlinear, width = 8, height = 6, dpi = 300)

# HS: 使用 rcs() 函数对饮酒量进行限制三次样条变换，节点数设为3（可根据需要调整）
HS_model <- cph(
  Surv(follow_up_time, Hemorrhagic_stroke) ~ rcs(alcohol_pw, 3) + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate ,
  data = data
)
# 查看模型摘要
summary(HS_model)  # 输出模型的详细结果，包括各变量的系数和显著性
# 检查非线性项是否显著
anova(HS_model)  # 对模型进行方差分析，检查饮酒量的非线性项是否显著

# 绘制限制三次样条图
# 使用 Predict() 函数生成饮酒量与疾病风险的预测值
pred_HS <- Predict(HS_model, alcohol_pw)  # 生成预测值
# 将 Log HR 转换为 HR
pred_HS$yhat <- exp(pred_HS$yhat)  # 对 yhat 进行指数变换
pred_HS$lower <- exp(pred_HS$lower)  # 对置信区间下限进行指数变换
pred_HS$upper <- exp(pred_HS$upper)  # 对置信区间上限进行指数变换
# 找到最小风险剂量（饮酒量对应的HR最低点）
min_risk_index <- which.min(pred_HS$yhat)
min_risk_alcohol <- pred_HS$alcohol_pw[min_risk_index]
min_risk_hr <- pred_HS$yhat[min_risk_index]
# 找到风险显著增加的阈值（HR置信区间下限超过1的最小饮酒量）
threshold_index <- which(pred_HS$lower > 1)
if (length(threshold_index) > 0) {
  threshold_alcohol <- pred_HS$alcohol_pw[threshold_index[1]]
  threshold_hr <- pred_HS$yhat[threshold_index[1]]
} else {
  threshold_alcohol <- NA
  threshold_hr <- NA
}
# 输出结果
cat("最小风险剂量（饮酒量）:", min_risk_alcohol, "\n")
cat("对应的危险比（HR）:", min_risk_hr, "\n")
if (!is.na(threshold_alcohol)) {
  cat("风险显著增加的阈值（饮酒量）:", threshold_alcohol, "\n")
  cat("对应的危险比（HR）:", threshold_hr, "\n")
} else {
  cat("未找到风险显著增加的阈值（置信区间未超过1）。\n")
}
# 绘制限制三次样条图
HS_nonlinear <- ggplot(pred_HS, aes(x = alcohol_pw_std, y = yhat)) +
  geom_line() +  # 绘制曲线
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  # 添加置信区间
  geom_vline(xintercept = min_risk_alcohol, linetype = "dashed", color = "blue") +  # 最小风险剂量
  geom_vline(xintercept = threshold_alcohol, linetype = "dashed", color = "red") +  # 风险显著增加的阈值
  labs(
    title = "The relationship between alcohol intake and incident Hemorrhagic stroke risk.",
    x = "Alcohol intake",
    y = "HR"
  ) +
  theme_minimal()  # 使用简洁的主题
# 显示图形
print(HS_nonlinear)
ggsave("non-linear HS.pdf", plot = HS_nonlinear, width = 8, height = 6, dpi = 300)

# IS: 使用 rcs() 函数对饮酒量进行限制三次样条变换，节点数设为3（可根据需要调整）
IS_model <- cph(
  Surv(follow_up_time, Ischemic_stroke) ~ rcs(alcohol_pw_std, 3) + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate ,
  data = data
)
# 查看模型摘要
summary(IS_model)  # 输出模型的详细结果，包括各变量的系数和显著性
# 检查非线性项是否显著
anova(IS_model)  # 对模型进行方差分析，检查饮酒量的非线性项是否显著

# 绘制限制三次样条图
# 使用 Predict() 函数生成饮酒量与疾病风险的预测值
pred_IS <- Predict(IS_model, alcohol_pw_std)  # 生成预测值
# 将 Log HR 转换为 HR
pred_IS$yhat <- exp(pred_IS$yhat)  # 对 yhat 进行指数变换
pred_IS$lower <- exp(pred_IS$lower)  # 对置信区间下限进行指数变换
pred_IS$upper <- exp(pred_IS$upper)  # 对置信区间上限进行指数变换
# 保存图形到本地文件
IS_nonlinear <- ggplot(pred_IS) +  # 使用ggplot2绘制图形
  geom_line() +  # 绘制曲线
  labs(
    title = "The relationship between alcohol intake and incident Ischemic stroke risk.",
    x = "Alcohol intake",
    y = "HR"
  ) +
  theme_minimal()  # 使用简洁的主题
ggsave("non-linear IS.pdf", plot = IS_nonlinear, width = 8, height = 6, dpi = 300)

# TIA: 使用 rcs() 函数对饮酒量进行限制三次样条变换，节点数设为3（可根据需要调整）
TIA_model <- cph(
  Surv(follow_up_time, TIA) ~ rcs(alcohol_pw_std, 3) + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate ,
  data = data
)
# 查看模型摘要
summary(TIA_model)  # 输出模型的详细结果，包括各变量的系数和显著性
# 检查非线性项是否显著
anova(TIA_model)  # 对模型进行方差分析，检查饮酒量的非线性项是否显著

# 绘制限制三次样条图
# 使用 Predict() 函数生成饮酒量与疾病风险的预测值
pred_TIA <- Predict(TIA_model, alcohol_pw_std)  # 生成预测值
# 将 Log HR 转换为 HR
pred_TIA$yhat <- exp(pred_TIA$yhat)  # 对 yhat 进行指数变换
pred_TIA$lower <- exp(pred_TIA$lower)  # 对置信区间下限进行指数变换
pred_TIA$upper <- exp(pred_TIA$upper)  # 对置信区间上限进行指数变换
# 保存图形到本地文件
TIA_nonlinear <- ggplot(pred_TIA) +  # 使用ggplot2绘制图形
  geom_line() +  # 绘制曲线
  labs(
    title = "The relationship between alcohol intake and incident TIA risk.",
    x = "Alcohol intake",
    y = "HR"
  ) +
  theme_minimal()  # 使用简洁的主题
ggsave("non-linear TIA.pdf", plot = TIA_nonlinear, width = 8, height = 6, dpi = 300)