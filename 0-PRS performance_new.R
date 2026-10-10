# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\alcohol投稿\\16-peerj major revision\\peerj major revision 新增材料")
getwd()

data <- read.csv("data with prs3_r.csv")
prs_p <- read.csv("prs performance.csv")

##### 排除非白人人群、缺乏基因数据的人群
# 统计非白人的人数
no_white <- sum(data$ethnicity == 2)
cat("非白人的人数：", no_white, "\n")
data <- data[data$ethnicity != 2, ]
# 统计无PRS的人群
no_PRS <- sum(is.na(data$PRS) | data$PRS == "", na.rm = TRUE)
cat("无PRS的人数为：", no_PRS, "\n")
data <- data[!is.na(data$PRS) & data$PRS != "", ]


library(dplyr)
prs_p <- prs_p %>%
  filter(eid %in% data$eid)
# 将data中的PRS列、alcohol_pw列按照eid合并到prs_p中
prs_p_updated <- prs_p %>%
  left_join(data %>% select(eid, PRS, alcohol_pw), by = "eid", suffix = c("", "_data"))


##### PRS与酒精摄入量的线性关联性分析
# 1. 数据预处理：将字符型变量转化为因子
prs_p_updated$p31 <- as.factor(prs_p_updated$p31)          # 性别
prs_p_updated$p22000 <- as.factor(prs_p_updated$p22000)    # 基因芯片批次
prs_p_updated$p54_i0 <- as.factor(prs_p_updated$p54_i0)    # 评估中心
# 2. 定义协变量公式 (前10个PC + 年龄 + 性别 + 批次 + 中心)
covariates <- "p21022 + p31 + p22000 + p54_i0 + 
               p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + 
               p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10"
# 3. 构建基础模型 (Base Model: 仅含协变量)
base_formula <- as.formula(paste("alcohol_pw ~", covariates))
model_base <- lm(base_formula, data = prs_p_updated)
# 4. 构建全模型 (Full Model: 协变量 + PRS)
full_formula <- as.formula(paste("alcohol_pw ~ PRS +", covariates))
model_full <- lm(full_formula, data = prs_p_updated)
# 5. 提取核心指标
summary_full <- summary(model_full)
# 计算 Incremental R-squared
r2_base <- summary(model_base)$r.squared
r2_full <- summary_full$r.squared
incremental_r2 <- r2_full - r2_base
# 提取 PRS 的 Beta, SE 和 P-value
prs_stats <- summary_full$coefficients["PRS", ]
# 6. 打印结果
cat("--- PRS Performance Analysis ---\n")
cat("Incremental R2 (PRS contribution):", round(incremental_r2 * 100, 4), "%\n")
cat("PRS Beta:", round(prs_stats["Estimate"], 4), "\n")
cat("PRS P-value:", format.pval(prs_stats["Pr(>|t|)"], digits = 4), "\n")


##### 分位数分析
library(ggplot2)
library(dplyr)
# 1. 确保使用的是你的数据框 prs_p_updated
# 按 PRS 分成 10 组 (已经运行过这一行的话可以跳过)
prs_p_updated$prs_decile <- ntile(prs_p_updated$PRS, 10)
# 2. 计算每一组的平均酒精摄入量和标准误 (SE)
# 注意：这里把 df 换成了 prs_p_updated
decile_summary <- prs_p_updated %>%
  group_by(prs_decile) %>%
  summarise(
    mean_intake = mean(alcohol_pw, na.rm = TRUE),
    se_intake = sd(alcohol_pw, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'  # 保持结果干净
  )
# 3. 检查一下计算结果
print(decile_summary)
# 4. 绘制柱状图
ggplot(decile_summary, aes(x = factor(prs_decile), y = mean_intake)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_intake - se_intake, ymax = mean_intake + se_intake), width = 0.2) +
  labs(
    title = "Alcohol Consumption by PRS Deciles",
    subtitle = paste("P-value < 2.2e-16, Incremental R2 = 0.1335%"),
    x = "PRS Deciles (1 = Lowest, 10 = Highest Risk)",
    y = "Mean Alcohol Intake (Units/Week)"
  ) +
  theme_minimal()
# 假设你已经运行了之前的绘图代码并得到了 ggplot 对象
# 如果没有赋值，可以先执行：
p <- ggplot(decile_summary, aes(x = factor(prs_decile), y = mean_intake)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "white", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_intake - se_intake, ymax = mean_intake + se_intake), width = 0.2) +
  labs(
    title = "Alcohol Consumption by PRS Deciles",
    subtitle = "P-value < 2.2e-16, Incremental R2 = 0.1335%",
    x = "PRS Deciles (1 = Lowest, 10 = Highest Risk)",
    y = "Mean Alcohol Intake (Units/Week)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
# 保存为 PDF
# width 和 height 的单位默认是英寸 (Inches)，7x5 是论文插图的常用比例
ggsave("PRS_Performance_Deciles.pdf", plot = p, width = 7, height = 5, device = "pdf")


##### logistic 回归
# 1. 对 PRS 进行标准化
prs_p_updated$PRS_z <- scale(prs_p_updated$PRS)

# 2. 准备二分类结局变量 (以酒精依赖为例，假设变量名为 AUD_status, 1为病例, 0为对照)
# 如果你没有诊断数据，可以按饮酒量前10%定义为 Heavy_Drinker
prs_p_updated$Heavy_Drinker <- ifelse(prs_p_updated$alcohol_pw > quantile(prs_p_updated$alcohol_pw, 0.9, na.rm=T), 1, 0)

# 3. 运行 Logistic 回归
# 协变量保持一致：Age, Sex, PC1-10, Batch, Center
logit_model <- glm(Heavy_Drinker ~ PRS_z + p21022 + p31 + p22000 + p54_i0 + 
                     p22009_a1 + p22009_a2 + p22009_a3 + p22009_a4 + p22009_a5 + 
                     p22009_a6 + p22009_a7 + p22009_a8 + p22009_a9 + p22009_a10, 
                   data = prs_p_updated, family = binomial)

# 4. 提取 OR 值和 95% 置信区间
summary_logit <- summary(logit_model)
OR <- exp(coef(logit_model)["PRS_z"])
CI <- exp(confint.default(logit_model)["PRS_z", ])
P_val <- summary_logit$coefficients["PRS_z", "Pr(>|z|)"]

# 5. 打印结果
cat("OR per SD of PRS:", round(OR, 3), "\n")
cat("95% CI:", round(CI[1], 3), "-", round(CI[2], 3), "\n")
cat("P-value:", format.pval(P_val, digits = 4), "\n")