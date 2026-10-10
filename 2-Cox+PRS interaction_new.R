# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\5-Cox")
getwd()
data <- read.csv("data with prs3_r.csv")

# 加载必要的R包
library(survival)    # Cox回归分析
library(caret)
library(broom)

# 变量类别处理
data$ACDs[data$ACDs == 2] <- 0
data$All_cause_stroke[is.na(data$All_cause_stroke)] <- 0
data$Hemorrhagic_stroke[is.na(data$Hemorrhagic_stroke)] <- 0
data$Ischemic_stroke[is.na(data$Ischemic_stroke)] <- 0
data$TIA_[is.na(data$TIA)] <- 0
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
data$alcohol_pw_std <- scale(data$alcohol_pw)  # 标准化处理
data$alcohol_pw_std <as.numeric(data$alcohol_pw_std)
str(data)

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


##### 排除非白人人群、缺乏基因数据的人群
# 统计非白人的人数
no_white <- sum(data$ethnicity == 2)
cat("非白人的人数：", no_white, "\n")
data <- data[data$ethnicity != 2, ]
# 统计无PRS的人群
no_PRS <- sum(is.na(data$PRS) | data$PRS == "", na.rm = TRUE)
cat("无PRS的人数为：", no_PRS, "\n")
data <- data[!is.na(data$PRS) & data$PRS != "", ]

##### 根据PRS评分的三分位数，将PRS风险分为3个等级（1-Low tertile, 2-Middle tertile, 3-High tertile）
## 直接将PRS列进行分析，得出的95%CI非常大
## 因此对PRS列进行处理：对数转换（log(1 + PRS)），后再将PRS列标准化（减去均值并除以标准差），使其均值为0，标准差为1
summary(data$PRS)
data$PRS_log <- log(1 + data$PRS)  # 对数转换
summary(data$PRS_log)
data$PRS_scaled <- scale(data$PRS_log)  # 标准化处理
summary(data$PRS_scaled)
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

quantiles <- quantile(data$PRS_scaled, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
cat("PRS列的三分位数为：\n")
print(quantiles)
# 新建PRS_risk，根据三分位数划分基因风险人群
# 计算PRS列的25%分位数和75%分位数
quantiles <- quantile(data$PRS_scaled, probs = c(0.25, 0.75), na.rm = TRUE)
Q1 <- quantiles[1]  # 25%分位数
Q3 <- quantiles[2]  # 75%分位数
cat("PRS列的25%分位数（Q1）：", Q1, "\n")
cat("PRS列的75%分位数（Q3）：", Q3, "\n")
data$PRS_risk <- ifelse(data$PRS_scaled <= Q1, "1",
                        ifelse(data$PRS_scaled > Q3, "3", "2"))
data$PRS_risk <- as.factor(data$PRS_risk)



##### Step 1 独立效应分析
library(survival)
# 定义生存分析的结果变量：假设你的数据中有 follow_up_time (生存时间) 和 ACDs (结局: 1=ACDs, 0=Censored)
# 模型 1 的协变量
cov_m1 <- c("preference", "with_food")
# 模型 2 的协变量 (在 M1 基础上增加人口学和生活方式)
cov_m2 <- c(cov_m1, "age", "sex", "ethnicity", "education", "TDI", "smoking", "diet_quality", "sleep_quality_mer")
# 模型 3 的协变量 (在 M2 基础上增加代谢与健康指标)
cov_m3 <- c(cov_m2, "blood_pressure_mer", "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")
# --- Model 1: 基础调整 ---
formula_m1 <- as.formula(paste("Surv(follow_up_time, ACDs) ~ alcohol_risk + PRS_risk +", paste(cov_m1, collapse = "+")))
fit_m1 <- coxph(formula_m1, data = data)
summary(fit_m1)
# --- Model 2: 中度调整 ---
formula_m2 <- as.formula(paste("Surv(follow_up_time, ACDs) ~ alcohol_risk + PRS_risk +", paste(cov_m2, collapse = "+")))
fit_m2 <- coxph(formula_m2, data = data)
summary(fit_m2)
# --- Model 3: 深度调整 (全模型) ---
formula_m3 <- as.formula(paste("Surv(follow_up_time, ACDs) ~ alcohol_risk + PRS_risk +", paste(cov_m3, collapse = "+")))
fit_m3 <- coxph(formula_m3, data = data)
summary(fit_m3)
# 提取结果的函数
extract_res <- function(model_fit, model_name) {
  res <- summary(model_fit)$conf.int
  p_val <- summary(model_fit)$coefficients
  # 提取 alcohol 和 PRS 的行
  rows <- grep("alcohol_risk|PRS_risk", rownames(res))
  data.frame(
    Model = model_name,
    Variable = rownames(res)[rows],
    HR = res[rows, 1],
    Lower_CI = res[rows, 3],
    Upper_CI = res[rows, 4],
    P_value = p_val[rows, 5]
  )
}
# 合并结果
step1_results <- rbind(
  extract_res(fit_m1, "Model 1"),
  extract_res(fit_m2, "Model 2"),
  extract_res(fit_m3, "Model 3")
)
print(step1_results)



##### Step 2 交互项 Cox 模型
library(survival)
# 1. 重新定义全模型协变量（延续 Model 3）
cov_m3 <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
            "TDI", "smoking", "diet_quality", "sleep_quality_mer", "blood_pressure_mer", 
            "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")
# 2. 构建交互作用公式
# 我们考察 alcohol_pw_std 与 PRS_scaled (连续型) 的交互
formula_int <- as.formula(paste(
  "Surv(follow_up_time, ACDs) ~ alcohol_risk * PRS_risk +", 
  paste(cov_m3, collapse = "+")
))
# 3. 拟合交互模型
fit_interation <- coxph(formula_int, data = data)
# 4. 查看结果
summary_int <- summary(fit_interation)
print(summary_int)
# 提取交互项的系数、HR 和 P 值
int_coefs <- summary_int$coefficients
int_rows <- grep(":", rownames(int_coefs)) # 匹配包含冒号的行
interaction_tab <- data.frame(
  Term = rownames(int_coefs)[int_rows],
  HR = exp(int_coefs[int_rows, "coef"]),
  Lower_CI = exp(int_coefs[int_rows, "coef"] - 1.96 * int_coefs[int_rows, "se(coef)"]),
  Upper_CI = exp(int_coefs[int_rows, "coef"] + 1.96 * int_coefs[int_rows, "se(coef)"]),
  P_value = int_coefs[int_rows, "Pr(>|z|)"]
)
print(interaction_tab)



##### Step 3 计算PRS与alcohol intake的相关性
# 不含 PRS
lm_base <- lm(alcohol_pw_std ~ + preference + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = data)
# 含 PRS
lm_full <- lm(alcohol_pw_std ~ PRS_scaled + preference + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = data)
# R²变化
r2_base <- summary(lm_base)$r.squared
r2_full <- summary(lm_full)$r.squared
delta_r2 <- r2_full - r2_base
delta_r2