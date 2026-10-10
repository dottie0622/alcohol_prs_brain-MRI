# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mri_mediation with time.csv")

# 加载必要的R包
library(ggplot2)
library(survival)    # Cox回归分析
library(lattice)
library(caret)
library(car)         # 计算VIF
library(carData)     # 计算VIF
library(broom)
library(dplyr)
library(tidyverse)
library(MASS)
library(Matrix)
library(mvtnorm)
library(sandwich)
library(mediation)
library(boot)


# 数据类型转换
mri <- mri %>%
  mutate(preference = case_when(
    preference == "wine" ~ "1",
    preference == "beer" ~ "2",
    preference == "spirits" ~ "3",
    preference == "mixed" ~ "4",
    TRUE ~ preference
  ))
mri <- mri %>%
  mutate(frequency = case_when(
    frequency == "Once or twice a week" ~ "1",
    frequency == "Three or four times a week" ~ "2",
    frequency == "Daily or almost daily" ~ "3",
    TRUE ~ frequency
  ))
mri <- mri %>%
  mutate(with_food = case_when(
    with_food == "Yes" ~ "1",
    with_food == "No" ~ "2",
    with_food == "It varies" ~ "3",
    TRUE ~ with_food
  ))
mri <- mri %>%
  mutate(ACDs = case_when(
    ACDs == 2 ~ 0,
    ACDs == 1 ~ 1,
    TRUE ~ ACDs
  ))
mri$alcohol_risk <- as.factor(mri$alcohol_risk)
mri$frequency <- as.factor(mri$frequency)
mri$preference <- as.factor(mri$preference)
mri$with_food <- as.factor(mri$with_food)
mri$sex <- as.factor(mri$sex)
mri$ethnicity <- as.factor(mri$ethnicity)
mri$education <- as.factor(mri$education)
mri$smoking <- as.factor(mri$smoking)
mri$diet_quality <- as.factor(mri$diet_quality)
mri$sleep_quality <- as.factor(mri$sleep_quality)
mri$blood_pressure <- as.factor(mri$blood_pressure)
mri$BMI <- as.factor(mri$BMI)
mri$health_rating <- as.factor(mri$health_rating)
mri$ACDs <- as.numeric(mri$ACDs)
mri$All_cause_stroke <- as.factor(mri$All_cause_stroke)
mri$Hemorrhagic_stroke <- as.factor(mri$Hemorrhagic_stroke)
mri$Ischemic_stroke <- as.factor(mri$Ischemic_stroke)
mri$TIA <- as.factor(mri$TIA)
str(mri)


##### 中介分析
# 步骤1: 拟合模型 中介变量 ~ 自变量
# 步骤2: 拟合模型 结局变量 ~ 自变量 + 中介变量
# 步骤3: 进行中介分析

# 自变量、协变量
mri$log_alcohol_pw <- log(mri$alcohol_pw + 1)  # 加1以避免对数为负无穷
exposure <- "log_alcohol_pw"
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality", "blood_pressure", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")

# 中介变量列表
mediators <- c("p25020_i2","p25016_i2","p25507_i2","p25502_i2","p25490_i2",
               "p25512_i2","p25508_i2","p25503_i2","p25499_i2","p25491_i2",
               "p25540_i2","p25538_i2","p25536_i2","p25531_i2","p25529_i2",
               "p25527_i2","p25523_i2","p25519_i2","p25517_i2","p25541_i2",
               "p25539_i2","p25537_i2","p25535_i2","p25530_i2","p25528_i2",
               "p25526_i2","p25524_i2","p25518_i2")

# 保存结果
results <- data.frame()


# 清理数据，避免缺失值问题
mri_clean <- na.omit(mri[, c("eid", exposure, "ACDs", "follow_up_time", covariates, mediators)])
mri_clean$log_alcohol_pw <- as.numeric(mri_clean$log_alcohol_pw)
# 把清洗后的数据框保存为mediation analysis data_MRI.csv
# 保存 CSV 文件
write.csv(mri_clean, "mediation analysis data_MRI.csv", row.names = FALSE)

##### ——————————————————————————————————————————————————————————————————————————
##### 直接使用清洗后的数据框进行数据分析
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mediation analysis data_MRI.csv")
exposure <- "log_alcohol_pw"
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality", "blood_pressure", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")
# 中介变量列表
mediators <- c("p25020_i2","p25016_i2","p25507_i2","p25502_i2","p25490_i2",
               "p25512_i2","p25508_i2","p25503_i2","p25499_i2","p25491_i2",
               "p25540_i2","p25538_i2","p25536_i2","p25531_i2","p25529_i2",
               "p25527_i2","p25523_i2","p25519_i2","p25517_i2","p25541_i2",
               "p25539_i2","p25537_i2","p25535_i2","p25530_i2","p25528_i2",
               "p25526_i2","p25524_i2","p25518_i2")
# 数据类型修正
mri$preference <- as.factor(mri$preference)
mri$with_food <- as.factor(mri$with_food)
mri$sex <- as.factor(mri$sex)
mri$ethnicity <- as.factor(mri$ethnicity)
mri$education <- as.factor(mri$education)
mri$smoking <- as.factor(mri$smoking)
mri$diet_quality <- as.factor(mri$diet_quality)
mri$sleep_quality <- as.factor(mri$sleep_quality)
mri$blood_pressure <- as.factor(mri$blood_pressure)
mri$BMI <- as.factor(mri$BMI)
mri$health_rating <- as.factor(mri$health_rating)
mri$ACDs <- as.numeric(mri$ACDs)
str(mri)

# 中介效应分析
for (mediator in mediators) {
  cat("处理中介变量:", mediator, "\n")
  
  med_formula <- as.formula(paste(mediator, "~", paste(c(exposure, covariates), collapse = " + ")))
  out_formula <- as.formula(paste("ACDs ~", paste(c(exposure, mediator, covariates), collapse = " + ")))
  
  med_model <- tryCatch(lm(med_formula, data = mri), error = function(e) NULL)
  out_model <- tryCatch(glm(out_formula, family = binomial(link = "logit"), data = mri), error = function(e) NULL)
  
  if (is.null(med_model) | is.null(out_model)) {
    cat("跳过变量", mediator, "：模型拟合失败\n")
    next
  }
  
  mediation_result <- tryCatch(
    mediate(med_model, out_model, treat = exposure, mediator = mediator,
            boot = TRUE, sims = 1000),
    error = function(e) {
      cat("中介分析失败:", mediator, "错误信息：", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(mediation_result)) next
  
  result_row <- data.frame(
    mediator = mediator,
    ACME_est = mediation_result$d0,
    ACME_p = mediation_result$d0.p,
    ADE_est = mediation_result$z0,
    ADE_p = mediation_result$z0.p,
    Total_effect = mediation_result$tau.coef,
    Total_p = mediation_result$tau.p,
    Prop_mediated = mediation_result$n0,
    sig_ACME = ifelse(mediation_result$d0.p < 0.05, "*", "")
  )
  
  results <- rbind(results, result_row)
}

print(results)

# 保存 CSV 文件
write.csv(results, "mediation_results_standardized_2.csv", row.names = FALSE)

# 安装并加载qvalue包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)
# 计算FDR校正后的P值
fdr_ACEM_p <- qvalue(results$ACME_p, pi0 = 0.05, method = "fdr")
fdr_ADE_p <- qvalue(results$ADE_p, pi0 = 0.05, method = "fdr")
fdr_Total_p <- qvalue(results$Total_p, pi0 = 0.05, method = "fdr")
print(fdr_ACEM_p)
print(fdr_ADE_p)
print(fdr_Total_p)


##### ————————————————————————————————————————————————————————————————————————
##### 以上单一变量的中介分析得出的结果并不理想，拟进行按照脑区功能维度进行分析的
##### 降维处理（Z-score + 平均），将这 14 个组合指标作为中介变量，进行同时建模的中介分析
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mediation analysis data_MRI.csv")
exposure <- "log_alcohol_pw"
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality", "blood_pressure", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")
# 数据类型修正
mri$preference <- as.factor(mri$preference)
mri$with_food <- as.factor(mri$with_food)
mri$sex <- as.factor(mri$sex)
mri$ethnicity <- as.factor(mri$ethnicity)
mri$education <- as.factor(mri$education)
mri$smoking <- as.factor(mri$smoking)
mri$diet_quality <- as.factor(mri$diet_quality)
mri$sleep_quality <- as.factor(mri$sleep_quality)
mri$blood_pressure <- as.factor(mri$blood_pressure)
mri$BMI <- as.factor(mri$BMI)
mri$health_rating <- as.factor(mri$health_rating)
mri$ACDs <- as.numeric(mri$ACDs)
str(mri)

# ---------------------
# 1. 创建14个维度组合变量（Z标准化后取平均）
# ---------------------
mri$group1_volume <- rowMeans(scale(mri[c("p25020_i2", "p25016_i2")]))
mri$group2_longitudinal_fasciculus <- rowMeans(scale(mri[c("p25502_i2", "p25503_i2", "p25529_i2", "p25530_i2")]))
mri$group3_ant_thalamic_radiation <- rowMeans(scale(mri[c("p25490_i2", "p25491_i2", "p25517_i2", "p25518_i2")]))
mri$group4_sup_thalamic_radiation <- rowMeans(scale(mri[c("p25512_i2", "p25538_i2", "p25539_i2")]))
mri$group5_post_thalamic_radiation <- rowMeans(scale(mri[c("p25507_i2", "p25508_i2", "p25535_i2")]))
mri$group6_forceps_major <- rowMeans(scale(mri[c("p25499_i2", "p25526_i2")]))
mri$group7_uncinate_fasciculus <- rowMeans(scale(mri[c("p25540_i2", "p25541_i2")]))
mri$group8_sup_long_fasciculus <- rowMeans(scale(mri[c("p25536_i2", "p25537_i2")]))
mri$group9_middle_cerebellar_peduncle <- scale(mri[["p25531_i2"]])
mri$group10_inferior_front_occipital <- rowMeans(scale(mri[c("p25527_i2", "p25528_i2")]))
mri$group11_corticospinal_tract <- rowMeans(scale(mri[c("p25523_i2", "p25524_i2")]))
mri$group12_cingulum <- scale(mri[["p25519_i2"]])

# ---------------------
# 2. 中介分析：遍历每个 group_* 中介变量
# ---------------------
library(mediation)
results <- data.frame()

mediator_vars <- c("group1_volume", "group2_longitudinal_fasciculus", "group3_ant_thalamic_radiation",
                   "group4_sup_thalamic_radiation", "group5_post_thalamic_radiation", "group6_forceps_major",
                   "group7_uncinate_fasciculus", "group8_sup_long_fasciculus", "group9_middle_cerebellar_peduncle",
                   "group10_inferior_front_occipital", "group11_corticospinal_tract", "group12_cingulum")
mri$group1_volume <- as.numeric(mri$group1_volume)
mri$group2_longitudinal_fasciculus <- as.numeric(mri$group2_longitudinal_fasciculus)
mri$group3_ant_thalamic_radiation <- as.numeric(mri$group3_ant_thalamic_radiation)
mri$group4_sup_thalamic_radiation <- as.numeric(mri$group4_sup_thalamic_radiation)
mri$group5_post_thalamic_radiation <- as.numeric(mri$group5_post_thalamic_radiation)
mri$group6_forceps_major <- as.numeric(mri$group6_forceps_major)
mri$group7_uncinate_fasciculus <- as.numeric(mri$group7_uncinate_fasciculus)
mri$group8_sup_long_fasciculus <- as.numeric(mri$group8_sup_long_fasciculus)
mri$group9_middle_cerebellar_peduncle <- as.numeric(mri$group9_middle_cerebellar_peduncle)
mri$group10_inferior_front_occipital <- as.numeric(mri$group10_inferior_front_occipital)
mri$group11_corticospinal_tract <- as.numeric(mri$group11_corticospinal_tract)
mri$group12_cingulum <- as.numeric(mri$group12_cingulum)
for (mediator in mediator_vars) {
  cat("处理中介维度:", mediator, "\n")
  
  # 模型公式
  med_formula <- as.formula(paste(mediator, "~", paste(c(exposure, covariates), collapse = " + ")))
  out_formula <- as.formula(paste("ACDs ~", paste(c(exposure, mediator, covariates), collapse = " + ")))
  
  # 拟合模型
  med_model <- lm(med_formula, data = mri)
  out_model <- glm(out_formula, family = binomial(link = "logit"), data = mri)
  
  # 中介分析
  mediation_result <- mediate(med_model, out_model, treat = exposure, mediator = mediator,
                              boot = TRUE, sims = 1000)
  
  # 提取结果
  result_row <- data.frame(
    mediator = mediator,
    ACME_est = mediation_result$d0,
    ACME_p = mediation_result$d0.p,
    ADE_est = mediation_result$z0,
    ADE_p = mediation_result$z0.p,
    Total_effect = mediation_result$tau.coef,
    Total_p = mediation_result$tau.p,
    Prop_mediated = mediation_result$n0,
    sig_ACME = ifelse(mediation_result$d0.p < 0.05, "*", "")
  )
  
  results <- rbind(results, result_row)
}

# 输出结果
print(results)

# 保存 CSV 文件
write.csv(results, "mediation_results_脑结构分组.csv", row.names = FALSE)

library(qvalue)
# 计算FDR校正后的P值
fdr_ACEM_p <- qvalue(results$ACME_p, pi0 = 0.05, method = "fdr")
fdr_ADE_p <- qvalue(results$ADE_p, pi0 = 0.05, method = "fdr")
fdr_Total_p <- qvalue(results$Total_p, pi0 = 0.05, method = "fdr")
print(fdr_ACEM_p)
print(fdr_ADE_p)
print(fdr_Total_p)



##### ——————————————————————————————————————————————————————————————————————————
##### 负值过多，重新分析
##### 中介分析
# 步骤1: 拟合模型 中介变量 ~ 自变量
# 步骤2: 拟合模型 结局变量 ~ 自变量 + 中介变量
# 步骤3: 进行中介分析
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mediation analysis data_MRI.csv")
exposure <- "log_alcohol_pw"
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality", "blood_pressure", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate")
# 数据类型修正
mri$preference <- as.factor(mri$preference)
mri$with_food <- as.factor(mri$with_food)
mri$sex <- as.factor(mri$sex)
mri$ethnicity <- as.factor(mri$ethnicity)
mri$education <- as.factor(mri$education)
mri$smoking <- as.factor(mri$smoking)
mri$diet_quality <- as.factor(mri$diet_quality)
mri$sleep_quality <- as.factor(mri$sleep_quality)
mri$blood_pressure <- as.factor(mri$blood_pressure)
mri$BMI <- as.factor(mri$BMI)
mri$health_rating <- as.factor(mri$health_rating)
mri$ACDs <- as.numeric(mri$ACDs)
str(mri)

# 中介变量列表
mediators <- c("p25020_i2","p25016_i2","p25507_i2","p25502_i2","p25490_i2",
               "p25512_i2","p25508_i2","p25503_i2","p25499_i2","p25491_i2",
               "p25540_i2","p25538_i2","p25536_i2","p25531_i2","p25529_i2",
               "p25527_i2","p25523_i2","p25519_i2","p25517_i2","p25541_i2",
               "p25539_i2","p25537_i2","p25535_i2","p25530_i2","p25528_i2",
               "p25526_i2","p25524_i2","p25518_i2")
## 数据标准化处理
# Volume 的 2 列
volume_columns <- c("p25020_i2", "p25016_i2")
# FA 的 8 列
fa_columns <- c("p25507_i2", "p25502_i2", "p25490_i2", "p25512_i2", "p25508_i2", "p25503_i2", "p25499_i2", "p25491_i2")
# MD 的 28 列
md_columns <- c("p25540_i2", "p25538_i2", "p25536_i2", "p25531_i2", "p25529_i2", "p25527_i2", "p25523_i2", "p25519_i2", 
                "p25517_i2", "p25541_i2", "p25539_i2", "p25537_i2", "p25535_i2", "p25530_i2", "p25528_i2", "p25526_i2", 
                "p25524_i2", "p25518_i2")

# 对 Volume 的 2 列进行标准化
mri[volume_columns] <- scale(mri[volume_columns])
# 对 FA 的 8 列进行标准化
mri[fa_columns] <- scale(mri[fa_columns])
# 对 MD 的 28 列进行标准化
mri[md_columns] <- scale(mri[md_columns])
# 确认生存时间有效
summary(mri$follow_up_time)
# 确认事件变量编码正确
table(mri$ACDs)
# 确保每个协变量都是合适类型
str(mri[, c(exposure, mediators, covariates, "ACDs", "follow_up_time")])

# 保存结果
results <- data.frame()


for (mediator in mediators) {
  cat("处理中介变量:", mediator, "\n")
  
  # 构造公式
  med_formula <- as.formula(paste(mediator, "~", paste(c(exposure, covariates), collapse = " + ")))
  out_formula <- as.formula(paste("ACDs ~", paste(c(exposure, mediator, covariates), collapse = " + ")))
  
  # 拟合模型
  med_model <- tryCatch(lm(med_formula, data = mri), error = function(e) NULL)
  out_model <- tryCatch(glm(out_formula, family = binomial(link = "logit"), data = mri), error = function(e) NULL)
  
  if (is.null(med_model) | is.null(out_model)) {
    cat("跳过变量", mediator, "：模型拟合失败\n")
    next
  }
  
  # 中介分析
  mediation_result <- tryCatch(
    mediate(med_model, out_model, treat = exposure, mediator = mediator, boot = TRUE, sims = 1000),
    error = function(e) {
      cat("中介分析失败:", mediator, "错误信息：", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(mediation_result)) next
  
  result_row <- data.frame(
    mediator = mediator,
    ACME_est = mediation_result$d0,
    ACME_p = mediation_result$d0.p,
    ADE_est = mediation_result$z0,
    ADE_p = mediation_result$z0.p,
    Total_effect = mediation_result$tau.coef,
    Total_p = mediation_result$tau.p,
    Prop_mediated = mediation_result$n0
  )
  
  results <- rbind(results, result_row)
}

print(results)
write.csv(results, "mediation_results_glm_scaled.csv", row.names = FALSE)


# FDR校正
library(qvalue)
# 计算FDR校正后的P值
fdr_ACEM_p <- qvalue(results$ACME_p, pi0 = 0.05, method = "fdr")
fdr_ADE_p <- qvalue(results$ADE_p, pi0 = 0.05, method = "fdr")
fdr_Total_p <- qvalue(results$Total_p, pi0 = 0.05, method = "fdr")
print(fdr_ACEM_p)
print(fdr_ADE_p)
print(fdr_Total_p)


##### ——————————————————————————————————————————————————————————————————————————
# 单个变量的中介路径图
library(grid)
library(gridExtra)
library(readr)
library(dplyr)

# 读取中介分析结果数据
df <- read_csv("mediation_results_glm_scaled.csv")

# 筛选出显著中介变量（例如 ACME_p_adjusted < 0.05）
df_sig <- df %>% filter(ACME_p_adjusted < 0.05)

# 添加显著性星号函数
star_marker <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

# 输出路径图为 PDF
pdf("mediation_path_diagrams.pdf", width = 6, height = 4)

for (i in 1:nrow(df_sig)) {
  row <- df_sig[i, ]
  
  ## 路径标签
  beta_a <- sprintf("βₐ = %.4f%s", row$ACME_est, star_marker(row$ACME_p_adjusted))
  beta_b <- sprintf("β_b = %.4f%s", row$ADE_est, star_marker(row$ADE_p_adjusted))
  beta_c <- sprintf("β_DE = %.4f%s", row$Total_effect, star_marker(row$Total_p_adjusted))
  prop_label <- sprintf("Prop. = %.2f%%%s", row$Prop_mediated * 100, star_marker(row$ACME_p_adjusted))
  
  ## 变量标签
  label_exposure <- "Alcohol Intake"
  label_mediator <- row$mediator
  label_outcome <- "ACDs"
  
  ## 绘图
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 5)))
  
  text_box <- function(label, row, col, gp = gpar(fontsize = 12), box_color = "lightgray") {
    grid.rect(gp = gpar(fill = box_color, col = NA), 
              vp = viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(label, vp = viewport(layout.pos.row = row, layout.pos.col = col), gp = gp)
  }
  
  arrow_line <- function(x0, y0, x1, y1, label = "", pos = 0.5, gp = gpar()) {
    grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
               arrow = arrow(length = unit(0.15, "inches")), gp = gp)
    if (label != "") {
      grid.text(label, x = unit(x0 + (x1 - x0) * pos, "npc"),
                y = unit(y0 + (y1 - y0) * pos + 0.04, "npc"), gp = gpar(fontsize = 10))
    }
  }
  
  # 画节点
  text_box(label_exposure, 3, 1, box_color = "#FFE0E0")
  text_box(label_mediator, 2, 3, box_color = "#D0F0F7")
  text_box(label_outcome, 3, 5, box_color = "#E8EAF6")
  
  # 画路径
  arrow_line(0.22, 0.5, 0.5, 0.7, label = beta_a)
  arrow_line(0.5, 0.7, 0.78, 0.5, label = beta_b)
  arrow_line(0.22, 0.5, 0.78, 0.5, label = beta_c)
  
  # 中介比例
  grid.text(prop_label, x = unit(0.5, "npc"), y = unit(0.35, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
}

dev.off()


### 多个变量的中介路径图
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# 改为向量化函数
star_marker <- Vectorize(function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
})



# 读取数据
df <- read_csv("mediation_results_glm_scaled.csv")

# 选择显著中介变量（你可按需调整 n）
n_mediators <- 3
df_sig <- df %>%
  filter(ACME_p_adjusted < 0.05,
         Prop_mediated >= 0, Prop_mediated <= 1) %>%
  arrange(desc(abs(ACME_est))) %>%
  head(n_mediators)

# 设置节点信息
exposure <- "Alcohol\nIntake"
outcome <- "ACDs"
mediators <- df_sig$mediator

# 节点列表
nodes <- data.frame(
  name = c(exposure, mediators, outcome),
  type = c("exposure", rep("mediator", length(mediators)), "outcome")
)

# 边列表（路径）
edges <- data.frame(
  from = c(rep(exposure, length(mediators)), mediators),
  to = c(mediators, rep(outcome, length(mediators))),
  label = c(
    paste0("ACME: ", sprintf("%.4f", df_sig$ACME_est), star_marker(df_sig$ACME_p_adjusted)),
    paste0("ADE: ", sprintf("%.4f", df_sig$ADE_est), star_marker(df_sig$ADE_p_adjusted))
  )
)

# 坐标布局
nodes$y <- ifelse(nodes$type == "exposure", 0,
                  ifelse(nodes$type == "outcome", 0, seq(1, -1, length.out = length(mediators))))
nodes$x <- ifelse(nodes$type == "exposure", 0,
                  ifelse(nodes$type == "outcome", 2, 1))

# 合并边信息
edges <- edges %>%
  left_join(nodes, by = c("from" = "name")) %>%
  rename(x_start = x, y_start = y) %>%
  left_join(nodes, by = c("to" = "name")) %>%
  rename(x_end = x, y_end = y)

# 绘图
p <- ggplot() +
  # 箭头线
  geom_segment(data = edges, aes(x = x_start, y = y_start,
                                 xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.25, "cm")),
               size = 0.6, color = "gray30") +
  # 路径文字
  geom_text(data = edges, aes(x = (x_start + x_end)/2,
                              y = (y_start + y_end)/2 + 0.1,
                              label = label), size = 3.2) +
  # 节点
  geom_label(data = nodes, aes(x = x, y = y, label = name,
                               fill = type), size = 4, label.size = 0.4) +
  scale_fill_manual(values = c(exposure = "#FBB4AE", mediator = "#B3CDE3", outcome = "#CCEBC5")) +
  theme_void() +
  theme(legend.position = "none")
print(p)
# 使用 ggsave 保存为 PDF 文件
ggsave("mediation analysis intergrated.pdf", plot = p, device = "pdf", width = 8, height = 6)