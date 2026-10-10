# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\5-Cox")
getwd()
data <- read.csv("data with prs3_r.csv")

# 加载必要的R包
library(survival)    # Cox回归分析
library(caret)
library(car)         # 计算VIF
library(carData)     # 计算VIF
library(broom)

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
data$alcohol_risk <- relevel(data$alcohol_risk, ref = "1")

# 检查酒精摄入与模型 2 中 10 个协变量对ACDs风险的相互作用
# 选择结果显著的协变量，单独进行的后续亚组分析

install.packages("data.table")
library(data.table)
data <- as.data.table(data)

##### alcohol_pw连续型变量与模型 2 中 10 个协变量的交互结果
# 创建一个空的数据框来存储结果
interaction_pw <- data.frame(协变量 = character(), P值 = numeric(), 符号 = character(), stringsAsFactors = FALSE)

# 列出所有协变量
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", 
                 "education", "TDI", "smoking", "diet_quality", "sleep_quality")
 
# 遍历每个协变量，计算与alcohol_pw的交互作用的P值和符号
for (covariate in covariates) {
   # 构建包含交互项的模型公式
   formula <- as.formula(paste("ACDs ~ alcohol_pw *", covariate))
   
   # 拟合模型
   model <- lm(formula, data = data)
   
   # 提取模型摘要
   model_summary <- summary(model)
   
   # 找到交互项的名称
   coef_names <- rownames(model_summary$coefficients)
   interaction_terms <- grep(paste0("alcohol_pw:", covariate), coef_names, value = TRUE)
   
   if (length(interaction_terms) > 0) {
     # 提取交互项的P值和系数
     interaction_p_values <- sapply(interaction_terms, function(term) {
       model_summary$coefficients[term, "Pr(>|t|)"]
     })
     interaction_coefs <- sapply(interaction_terms, function(term) {
       coef(model)[term]
     })
     
     # 计算整体P值（取最小P值）
     overall_p_value <- min(interaction_p_values)
     
     # 根据系数的符号分配符号（取出现次数最多的符号）
     signs <- ifelse(interaction_coefs > 0, "+", ifelse(interaction_coefs < 0, "-", "0"))
     overall_sign <- names(which.max(table(signs)))
     
     # 将结果添加到数据框中
     interaction_pw <- rbind(interaction_pw, data.frame(协变量 = covariate, P值 = overall_p_value, 符号 = overall_sign))
   } else {
     # 如果没有找到交互项，添加一个NA行
     interaction_pw <- rbind(interaction_pw, data.frame(协变量 = covariate, P值 = NA, 符号 = NA))
   }
 }

# 查看结果
print(interaction_pw)
write.csv(interaction_pw, "alc-pw interaction p-value model2.csv", row.names = FALSE, quote = FALSE )

##### alcohol_risk分类变量与模型 2 中 10 个协变量的交互结果
# 加载必要的库
library(dplyr)

# 创建一个空的数据框来存储结果
interaction_risk <- data.frame(协变量 = character(), 交互水平 = character(), P值 = numeric(), 符号 = character(), stringsAsFactors = FALSE)

# 列出所有协变量
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", 
                "education", "TDI", "smoking", "diet_quality", "sleep_quality")

# 获取 alcohol_risk 的所有水平
levels_alcohol_risk <- levels(factor(data$alcohol_risk))

# 遍历每个协变量，计算与 alcohol_risk 的每个分类的交互作用的 P 值和符号
for (covariate in covariates) {
  for (level in 2:length(levels_alcohol_risk)) {  # 只考虑非参考水平，即2和3
    # 构建包含交互项的模型公式
    formula <- as.formula(paste("ACDs ~ alcohol_risk *", covariate))
    
    # 拟合模型
    model <- lm(formula, data = data)
    
    # 提取模型摘要
    model_summary <- summary(model)
    
    # 获取所有交互项系数的名称
    coef_names <- rownames(model_summary$coefficients)
    
    # 避免多次列出同一个交互项，匹配 alcohol_risk 与协变量的交互项
    interaction_terms <- grep(paste0("alcohol_risk", level, ":"), coef_names, value = TRUE)
    
    # 只保留主要交互项（即只保留 alcohol_risk 与协变量的交互项，不考虑因子水平）
    interaction_terms <- interaction_terms[grepl(covariate, interaction_terms)]
    
    # 如果模型中存在交互项
    if (length(interaction_terms) > 0) {
      # 提取交互项的 P 值和系数
      interaction_p_values <- sapply(interaction_terms, function(term) {
        model_summary$coefficients[term, "Pr(>|t|)"]
      })
      interaction_coefs <- sapply(interaction_terms, function(term) {
        coef(model)[term]
      })
      
      # 找到具有最小 P 值的交互项
      min_p_index <- which.min(interaction_p_values)
      interaction_p_value <- interaction_p_values[min_p_index]
      interaction_coef <- interaction_coefs[min_p_index]
      
      # 根据系数的符号分配符号
      if (interaction_coef > 0) {
        sign <- "+"
      } else if (interaction_coef < 0) {
        sign <- "-"
      } else {
        sign <- "0"
      }
      
      # 将结果添加到数据框中
      interaction_risk <- bind_rows(interaction_risk, data.frame(协变量 = covariate, 交互水平 = paste("alcohol_risk", level), P值 = interaction_p_value, 符号 = sign))
    } else {
      # 如果没有找到交互项，添加一个 NA 行
      interaction_risk <- bind_rows(interaction_risk, data.frame(协变量 = covariate, 交互水平 = paste("alcohol_risk", level), P值 = NA, 符号 = NA))
    }
  }
}

# 查看结果
print(interaction_risk)
write.csv(interaction_risk, "alc-risk interaction p-value of model2.csv", row.names = FALSE, quote = FALSE )


##### 根据上述交互分析得出了7个结果显著的协变量
##### 7个亚组协变量：preference/with_food/age/sex/education/diet_quality/sleep_quality
##### 选择其中4个所有P值均<0.05的协变量进行亚组分析：preference/with_food/age/sex
data <- data %>%
  mutate(sleep_quality_mer = case_when(
    sleep_quality %in% c(1, 2) ~ 1,
    sleep_quality == 3 ~ 2,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
data$sleep_quality_mer <- as.factor(data$sleep_quality_mer)
data <- data %>%
  mutate(blood_pressure_mer = case_when(
    blood_pressure %in% c(1, 2) ~ 1,
    blood_pressure == 3 ~ 2,
    blood_pressure == 4 ~ 3,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
data$blood_pressure_mer <- as.factor(data$blood_pressure_mer)

data$alcohol_pw_std <- scale(data$alcohol_pw)  

### preference亚组
summary(data$preference)
# 按preference分组
data_wine <- data %>% filter(preference == "1")
data_beer <- data %>% filter(preference == "2")
data_spirits <- data %>% filter(preference == "3")
data_mixed <- data %>% filter(preference == "4")
##### 在不同组别中的拟合 Cox 模型
#alcohol_pw_std与ACDs
cox_wine <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)
# alcohol_risk与ACDs
cox_wine <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)

#alcohol_pw_std与AS
cox_wine <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)
# alcohol_risk与ACDs
cox_wine <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)

# alcohol_pw_std与Hemorrhagic stroke
cox_wine <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)
# alcohol_risk与Hemorrhagic_stroke
cox_wine <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)

# alcohol_pw_std与Ischemic stroke
cox_wine <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)
# alcohol_risk与Ischemic stroke
cox_wine <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)

# alcohol_pw_std与TIA
cox_wine <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)
# alcohol_risk与TIA
cox_wine <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_wine)
summary(cox_wine)
cox_beer <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_beer)
summary(cox_beer)
cox_spirits <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_spirits)
summary(cox_spirits)
cox_mixed <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + with_food + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_mixed)
summary(cox_mixed)

### with_food亚组
summary(data$with_food)
# 按with_food分组
data_yes <- data %>% filter(with_food == "1")
data_no <- data %>% filter(with_food == "2")
data_varies <- data %>% filter(with_food == "3")
#alcohol_pw_std与ACDs
cox_yes <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)
# alcohol_risk与ACDs
cox_yes <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)

#alcohol_pw_std与AS
cox_yes <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)
# alcohol_risk与AS
cox_yes <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)

#alcohol_pw_std与Hemorrhagic_stroke
cox_yes <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)
# alcohol_risk与Hemorrhagic_stroke
cox_yes <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)

#alcohol_pw_std与Ischemic_stroke
cox_yes <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)
# alcohol_risk与Hemorrhagic_stroke
cox_yes <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)

#alcohol_pw_std与TIA
cox_yes <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)
# alcohol_risk与Hemorrhagic_stroke
cox_yes <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_yes)
summary(cox_yes)
cox_no <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_no)
summary(cox_no)
cox_varies <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + age + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_varies)
summary(cox_varies)

### 年龄亚组（＜65 years，≥65 years）
summary(data$age)
# 按age分组
data_elder <- data %>% filter(age > 64)
data_young <- data %>% filter(age < 65)
#alcohol_pw_std与ACDs
cox_young <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)
# alcohol_risk与ACDs
cox_young <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)

#alcohol_pw_std与AS
cox_young <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)
# alcohol_risk与AS
cox_young <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)

#alcohol_pw_std与Hemorrhagic_stroke
cox_young <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)
# alcohol_risk与Hemorrhagic_stroke
cox_young <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)

#alcohol_pw_std与Ischemic_stroke
cox_young <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)
# alcohol_risk与Ischemic_stroke
cox_young <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)

#alcohol_pw_std与TIA
cox_young <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)
# alcohol_risk与TIA
cox_young <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_young)
summary(cox_young)
cox_elder <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food + sex + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_elder)
summary(cox_elder)

### 性别亚组
summary(data$sex)
# 按sex分组
data_male <- data %>% filter(sex == "1")
data_female <- data %>% filter(sex == "2")
#alcohol_pw_std与ACDs
cox_male <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)
# alcohol_risk与ACDs
cox_male <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, ACDs) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)

#alcohol_pw_std与AS
cox_male <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)
# alcohol_risk与AS
cox_male <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, All_cause_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)

#alcohol_pw_std与Hemorrhagic_stroke
cox_male <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)
# alcohol_risk与Hemorrhagic_stroke
cox_male <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, Hemorrhagic_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)

#alcohol_pw_std与Ischemic_stroke
cox_male <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)
# alcohol_risk与Hemorrhagic_stroke
cox_male <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, Ischemic_stroke) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)

#alcohol_pw_std与TIA
cox_male <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, TIA) ~ alcohol_pw_std + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)
# alcohol_risk与TIA
cox_male <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_male)
summary(cox_male)
cox_female <- coxph(Surv(follow_up_time, TIA) ~ alcohol_risk + preference + with_food + age + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = data_female)
summary(cox_female)