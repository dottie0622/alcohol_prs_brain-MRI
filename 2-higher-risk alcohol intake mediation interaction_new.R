# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\alcohol投稿\\16-peerj major revision\\peerj major revision 新增材料")
getwd()

mri <- read.csv("mri_mediation with time.csv")


# 加载必要的R包
library(ggplot2)
library(survival)    # Cox回归分析
library(caret)
library(car)         # 计算VIF
library(carData)     # 计算VIF
library(broom)
library(dplyr)
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
str(mri)
mri$alcohol_pw <- scale(mri$alcohol_pw)
## 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
## 由于由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并（将1和2合并为1-poor or moderate，3改为2-healthy）
## blood_pressure3、blood_pressure4的VIF大于5 需要进行成分合并（将1和2合并为1-normal or elevated，3改为2-stage 1, 4改为3-stage 2）
table(mri$blood_pressure)
chisq.test(table(mri$blood_pressure))
mri <- mri %>%
  mutate(blood_pressure_mer = case_when(
    blood_pressure %in% c(1, 2) ~ 1,
    blood_pressure == 3 ~ 2,
    blood_pressure == 4 ~ 3,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
mri$blood_pressure_mer <- as.factor(mri$blood_pressure_mer)
## 由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并
table(mri$sleep_quality)
chisq.test(table(mri$sleep_quality))
mri <- mri %>%
  mutate(sleep_quality_mer = case_when(
    sleep_quality %in% c(1, 2) ~ 1,
    sleep_quality == 3 ~ 2,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
mri$sleep_quality_mer <- as.factor(mri$sleep_quality_mer)
mri$ACDs[mri$ACDs == 2] <- 0
sum(mri$Ischemic_stroke == 1, na.rm = TRUE)
##### 模型1
##### 皮层下区域
# 检查模型1协变量之间的多重共线性
model_1 <- vif(lm(alcohol_risk ~ preference + with_food, data = mri))
print(model_1)
# accumbens_left
model1_1L <- lm(p25023_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_1L)
# accumbens_right
model1_1R <- lm(p25024_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_1R)
# amygdala_left
model1_2L <- lm(p25021_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_2L)
# amygdala_right
model1_2R <- lm(p25022_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_2R)
# caudate_left
model1_3L <- lm(p25013_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_3L)
# caudate_right
model1_3R <- lm(p25014_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_3R)
# hippocampus_left
model1_4L <- lm(p25019_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_4L)
# hippocampus_right
model1_4R <- lm(p25020_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_4R)
# pallidum_left
model1_5L <- lm(p25017_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_5L)
# pallidum_right
model1_5R <- lm(p25018_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_5R)
# putamen_left
model1_6L <- lm(p25015_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_6L)
# putamen_right
model1_6R <- lm(p25016_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_6R)
# thalamus_left
model1_7L <- lm(p25011_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_7L)
# thalamus_right
model1_7R <- lm(p25012_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_7R)


##### 各个脑区的体积
# cuneus_left
model1_8L <- lm(p27207_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_8L)
# cuneus_right
model1_8R <- lm(p27300_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_8R)
# entorhinal_left
model1_9L <- lm(p27208_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_9L)
# entorhinal_right
model1_9R <- lm(p27301_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_9R)
# fusiform_left
model1_10L <- lm(p27209_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_10L)
# fusiform_right
model1_10R <- lm(p27302_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_10R)
# inferiorparietal_left
model1_11L <- lm(p27210_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_11L)
# inferiorparietal_right
model1_11R <- lm(p27303_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_11R)
# inferiortemporal_left
model1_12L <- lm(p27211_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_12L)
# inferiortemporal_right
model1_12R <- lm(p27304_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_12R)
# insula_left
model1_13L <- lm(p27235_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_13L)
# insula_right
model1_13R <- lm(p27328_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_13R)
# isthmuscingulate_left
model1_14L <- lm(p27212_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_14L)
# isthmuscingulate_right
model1_14R <- lm(p27305_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_14R)
# lateraloccipital_left
model1_15L <- lm(p27213_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_15L)
# lateraloccipital_right
model1_15R <- lm(p27306_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_15R)
# lateralorbitofrontal_left
model1_16L <- lm(p27214_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_16L)
# lateralorbitofrontal_right
model1_16R <- lm(p27307_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_16R)
# lingual_left
model1_17L <- lm(p27215_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_17L)
# lingual_right
model1_17R <- lm(p27308_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_17R)
# medialorbitofrontal_left
model1_18L <- lm(p27216_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_18L)
# medialorbitofrontal_right
model1_18R <- lm(p27309_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_18R)
# middletemporal_left
model1_19L <- lm(p27217_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_19L)
# middletemporal_right
model1_19R <- lm(p27310_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_19R)
# paracentral_left
model1_20L <- lm(p27219_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_20L)
# paracentral_right
model1_20R <- lm(p27312_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_20R)
# parahippocampal_left
model1_21L <- lm(p27218_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_21L)
# parahippocampal_right
model1_21R <- lm(p27311_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_21R)
# parsopercularis_left
model1_22L <- lm(p27220_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_22L)
# parsopercularis_right
model1_22R <- lm(p27313_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_22R)
# parsorbitalis_left
model1_23L <- lm(p27221_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_23L)
# parsorbitalis_right
model1_23R <- lm(p27314_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_23R)
# parstriangularis_left
model1_24L <- lm(p27222_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_24L)
# parstriangularis_right
model1_24R <- lm(p27315_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_24R)
# pericalcarine_left
model1_25L <- lm(p27223_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_25L)
# pericalcarine_right
model1_25R <- lm(p27316_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_25R)
# postcentral_left
model1_26L <- lm(p27224_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_26L)
# postcentral_right
model1_26R <- lm(p27317_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_26R)
# posteriorcingulate_left
model1_27L <- lm(p27225_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_27L)
# posteriorcingulate_right
model1_27R <- lm(p27318_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_27R)
# precentral_left
model1_28L <- lm(p27226_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_28L)
# precentral_right
model1_28R <- lm(p27319_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_28R)
# precuneus_left
model1_29L <- lm(p27227_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_29L)
# precuneus_right
model1_29R <- lm(p27320_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_29R)
# rostralanteriorcingulate_left
model1_30L <- lm(p27228_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_30L)
# rostralanteriorcingulate_right
model1_30R <- lm(p27321_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_30R)
# rostralmiddlefrontal_left
model1_31L <- lm(p27229_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_31L)
# rostralmiddlefrontal_right
model1_31R <- lm(p27322_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_31R)
# superiorfrontal_left
model1_32L <- lm(p27230_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_32L)
# superiorfrontal_right
model1_32R <- lm(p27323_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_32R)
# superiorparietal_left
model1_33L<- lm(p27231_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_33L)
# superiorparietal_right
model1_33R <- lm(p27324_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_33R)
# superiortemporal_left
model1_34L <- lm(p27232_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_34L)
# superiortemporal_right
model1_34R <- lm(p27325_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_34R)
# supramarginal_left
model1_35L <- lm(p27233_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_35L)
# supramarginal_right
model1_35R <- lm(p27326_i2 ~ alcohol_risk + preference + with_food, data = mri)
summary(model1_35R)
# transversetemporal_left
model1_36L <- lm(p27234_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_36L)
# transversetemporal_right
model1_36R <- lm(p27327_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_36R)
# caudalanteriorcingulate_left
model1_37L <- lm(p27205_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_37L)
# caudalanteriorcingulate_right
model1_37R <- lm(p27298_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_37R)
# caudalmiddlefrontal_left
model1_38L <- lm(p27206_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_38L)
# caudalmiddlefrontal_right
model1_38R <- lm(p27299_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer, data = mri)
summary(model1_38R)
# 提取每个模型的p值（以alcohol_risk为例）
p_values <- c(
  summary(model1_1L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # accumbens_left
  summary(model1_1R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # accumbens_right
  summary(model1_2L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # amygdala_left
  summary(model1_2R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # amygdala_right
  summary(model1_3L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # caudate_left
  summary(model1_3R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # caudate_right
  summary(model1_4L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # hippocampus_left
  summary(model1_4R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # hippocampus_right
  summary(model1_5L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # pallidum_left
  summary(model1_5R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # pallidum_right
  summary(model1_6L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # putamen_left
  summary(model1_6R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # putamen_right
  summary(model1_7L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # thalamus_left
  summary(model1_7R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # thalamus_right
  summary(model1_8L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # cuneus_left
  summary(model1_8R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # cuneus_right
  summary(model1_9L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # entorhinal_left
  summary(model1_9R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # entorhinal_right
  summary(model1_10L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # fusiform_left
  summary(model1_10R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # fusiform_right
  summary(model1_11L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # inferiorparietal_left
  summary(model1_11R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # inferiorparietal_right
  summary(model1_12L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # inferiortemporal_left
  summary(model1_12R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # inferiortemporal_right
  summary(model1_13L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # insula_left
  summary(model1_13R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # insula_right
  summary(model1_14L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # isthmuscingulate_left
  summary(model1_14R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # isthmuscingulate_right
  summary(model1_15L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lateraloccipital_left
  summary(model1_15R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lateraloccipital_right
  summary(model1_16L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lateralorbitofrontal_left
  summary(model1_16R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lateralorbitofrontal_right
  summary(model1_17L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lingual_left
  summary(model1_17R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # lingual_right
  summary(model1_18L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # medialorbitofrontal_left
  summary(model1_18R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # medialorbitofrontal_right
  summary(model1_19L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # middletemporal_left
  summary(model1_19R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # middletemporal_right
  summary(model1_20L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # paracentral_left
  summary(model1_20R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # paracentral_right
  summary(model1_21L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parahippocampal_left
  summary(model1_21R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parahippocampal_right
  summary(model1_22L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parsopercularis_left
  summary(model1_22R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parsopercularis_right
  summary(model1_23L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parsorbitalis_left
  summary(model1_23R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parsorbitalis_right
  summary(model1_24L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parstriangularis_left
  summary(model1_24R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # parstriangularis_right
  summary(model1_25L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # pericalcarine_left
  summary(model1_25R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # pericalcarine_right
  summary(model1_26L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # postcentral_left
  summary(model1_26R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # postcentral_right
  summary(model1_27L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # posteriorcingulate_left
  summary(model1_27R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # posteriorcingulate_right
  summary(model1_28L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # precentral_left
  summary(model1_28R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # precentral_right
  summary(model1_29L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # precuneus_left
  summary(model1_29R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # precuneus_right
  summary(model1_30L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # rostralanteriorcingulate_left
  summary(model1_30R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # rostralanteriorcingulate_right
  summary(model1_31L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # rostralmiddlefrontal_left
  summary(model1_31R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # rostralmiddlefrontal_right
  summary(model1_32L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiorfrontal_left
  summary(model1_32R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiorfrontal_right
  summary(model1_33L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiorparietal_left
  summary(model1_33R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiorparietal_right
  summary(model1_34L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiortemporal_left
  summary(model1_34R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # superiortemporal_right
  summary(model1_35L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # supramarginal_left
  summary(model1_35R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # supramarginal_right
  summary(model1_36L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # transversetemporal_left
  summary(model1_36R)$coefficients["alcohol_risk3", "Pr(>|t|)"], # transversetemporal_right
  summary(model1_37L)$coefficients["alcohol_risk3", "Pr(>|t|)"], # caudalanteriorcingulate_left
  summary(model1_37R)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # caudalanteriorcingulate_right
  summary(model1_38L)$coefficients["alcohol_risk3", "Pr(>|t|)"],  # caudalmiddlefrontal_left
  summary(model1_38R)$coefficients["alcohol_risk3", "Pr(>|t|)"]   # caudalmiddlefrontal_right
)
# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)




#####FA放大1000倍、MD放大10000倍，alcohol_pw标准化处理
names(mri)
mri[, 78:104] <- mri[, 78:104] * 1000
mri[, 105:131] <- mri[, 105:131] * 10000
##### 模型1
##### FA
# acoustic radiation_left
model1_1L <- lm(p25488_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_1L)
# acoustic radiation_right
model1_1R <- lm(p25489_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_1R)
# anterior thalamic radiation_left
model1_2L <- lm(p25490_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_2L)
# anterior thalamic radiation_right
model1_2R <- lm(p25491_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_2R)
# cingulate gyrus part of cingulum_left
model1_3L <- lm(p25492_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_3L)
# cingulate gyrus part of cingulum_right
model1_3R <- lm(p25493_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_3R)
# corticospinal tract_left
model1_4L <- lm(p25496_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_4L)
# corticospinal tract_right
model1_4R <- lm(p25497_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_4R)
# forceps major_left
model1_5L <- lm(p25498_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_5L)
# forceps major_right
model1_5R <- lm(p25499_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_5R)
# inferior fronto-occipital fasciculus_left
model1_6L <- lm(p25500_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_6L)
# inferior fronto-occipital fasciculus_right
model1_6R <- lm(p25501_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_6R)
# longitudinal fasciculus_left
model1_7L <- lm(p25502_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_7L)
# longitudinal fasciculus_right
model1_7R <- lm(p25503_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_7R)
# medial lemniscus_left
model1_8L <- lm(p25505_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_8L)
# medial lemniscus_right
model1_8R <- lm(p25506_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_8R)
# middle cerebellar peduncle
model1_9 <- lm(p25504_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_9)
# parahippocampal part of cingulum_left
model1_10L <- lm(p25494_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_10L)
# parahippocampal part of cingulum_right
model1_10R <- lm(p25495_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_10R)
# posterior thalamic radiation_left
model1_11L <- lm(p25507_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_11L)
# posterior thalamic radiation_right
model1_11R <- lm(p25508_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_11R)
# superior longitudinal fasciculus_left
model1_12L <- lm(p25509_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_12L)
# superior longitudinal fasciculus_right
model1_12R <- lm(p25510_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_12R)
# superior thalamic radiation_left
model1_13L <- lm(p25511_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_13L)
# superior thalamic radiation_right
model1_13R <- lm(p25512_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_13R)
# uncinate fasciculus_left
model1_14L <- lm(p25513_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_14L)
# uncinate fasciculus_right
model1_14R <- lm(p25514_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_14R)
##### 模型1
##### MD
# acoustic radiation_left
model1_15L <- lm(p25515_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_15L)
# acoustic radiation_right
model1_15R <- lm(p25516_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_15R)
# anterior thalamic radiation_left
model1_16L <- lm(p25517_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_16L)
# anterior thalamic radiation_right
model1_16R <- lm(p25518_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_16R)
# cingulate gyrus part of cingulum_left
model1_17L <- lm(p25519_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_17L)
# cingulate gyrus part of cingulum_right
model1_17R <- lm(p25520_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_17R)
# corticospinal tract_left
model1_18L <- lm(p25523_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_18L)
# corticospinal tract_right
model1_18R <- lm(p25524_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_18R)
# forceps major_left
model1_19L <- lm(p25525_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_19L)
# forceps major_right
model1_19R <- lm(p25526_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_19R)
# inferior fronto-occipital fasciculus_left
model1_20L <- lm(p25527_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_20L)
# inferior fronto-occipital fasciculus_right
model1_20R <- lm(p25528_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_20R)
# longitudinal fasciculus_left
model1_21L <- lm(p25529_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_21L)
# longitudinal fasciculus_right
model1_21R <- lm(p25530_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_21R)
# medial lemniscus_left
model1_22L <- lm(p25532_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_22L)
# medial lemniscus_right
model1_22R <- lm(p25533_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_22R)
# middle cerebellar peduncle
model1_23 <- lm(p25531_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_23)
# parahippocampal part of cingulum_left
model1_24L <- lm(p25521_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_24L)
# parahippocampal part of cingulum_right
model1_24R <- lm(p25522_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_24R)
# posterior thalamic radiation_left
model1_25L <- lm(p25534_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_25L)
# posterior thalamic radiation_right
model1_25R <- lm(p25535_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_25R)
# superior longitudinal fasciculus_left
model1_26L <- lm(p25536_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_26L)
# superior longitudinal fasciculus_right
model1_26R <- lm(p25537_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_26R)
# superior thalamic radiation_left
model1_27L <- lm(p25538_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_27L)
# superior thalamic radiation_right
model1_27R <- lm(p25539_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_27R)
# uncinate fasciculus_left
model1_28L <- lm(p25540_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_28L)
# uncinate fasciculus_right
model1_28R <- lm(p25541_i2 ~ alcohol_risk + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality_mer + blood_pressure_mer + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model1_28R)



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


##### 中介分析
# 步骤1: 拟合模型 中介变量 ~ 自变量
# 步骤2: 拟合模型 结局变量 ~ 自变量 + 中介变量
# 步骤3: 进行中介分析

# 自变量、协变量
mri$alcohol_risk <- mri$alcohol_risk  
mri$p25781_i2 <- log(mri$p25781_i2)
exposure <- "alcohol_risk"
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality_mer", "blood_pressure_mer", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate", "p25781_i2")

# 中介变量列表
mediators <- c("p25013_i2","p25014_i2","p25020_i2","p25016_i2",
               "p25490_i2","p25491_i2","p25499_i2","p25507_i2","p25508_i2","p25513_i2",
               "p25515_i2","p25516_i2","p25517_i2","p25518_i2",
               "p25523_i2","p25524_i2","p25526_i2","p25527_i2","p25528_i2","p25529_i2",
               "p25530_i2","p25532_i2","p25533_i2","p25531_i2","p25534_i2","p25535_i2",
               "p25536_i2","p25537_i2","p25538_i2","p25539_i2","p25540_i2","p25541_i2")

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
setwd("D:\\读研倒计时\\课题\\UKB\\data\\alcohol投稿\\16-peerj major revision\\peerj major revision 新增材料")
getwd()
mri <- read.csv("mri_mediation with time.csv")

# 1. 定义各组变量名
exposure <- "alcohol_risk"
outcome <- "ACDs"
potential_covs <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                    "TDI", "smoking", "diet_quality", "sleep_quality", "sleep_quality_mer",
                    "blood_pressure", "blood_pressure_mer", "BMI", "health_rating", 
                    "triglycerides", "HDL", "LDL", "glucose", "urate", "p25781_i2")
potential_meds <- c("p25013_i2","p25014_i2","p25020_i2", "p25016_i2",
                    "p25490_i2","p25491_i2","p25499_i2","p25507_i2","p25508_i2","p25513_i2",
                    "p25515_i2","p25516_i2","p25517_i2","p25518_i2",
                    "p25523_i2","p25524_i2","p25526_i2","p25527_i2","p25528_i2","p25529_i2",
                    "p25530_i2","p25532_i2","p25533_i2","p25531_i2","p25534_i2","p25535_i2",
                    "p25536_i2","p25537_i2","p25538_i2","p25539_i2","p25540_i2","p25541_i2")

# 2. 合并所有想要保留的列名
all_desired_columns <- c(exposure, outcome, potential_covs, potential_meds)
# 3. 核心步骤：只保留数据框中确实存在的列
# 使用 intersect 自动过滤掉数据框里没有的变量名，防止报错
existing_columns <- intersect(all_desired_columns, names(mri))
# 4. 执行筛选
mri_reduced <- mri[, existing_columns]
# 5. 验证结果
cat("原始列数：", ncol(mri), "\n")
cat("保留后的列数：", ncol(mri_reduced), "\n")
# 查看保留后的变量名
print(names(mri_reduced))


# 数据类型转换
mri_reduced <- mri_reduced %>%
  mutate(preference = case_when(
    preference == "wine" ~ "1",
    preference == "beer" ~ "2",
    preference == "spirits" ~ "3",
    preference == "mixed" ~ "4",
    TRUE ~ preference
  ))
mri_reduced <- mri_reduced %>%
  mutate(with_food = case_when(
    with_food == "Yes" ~ "1",
    with_food == "No" ~ "2",
    with_food == "It varies" ~ "3",
    TRUE ~ with_food
  ))
mri_reduced$alcohol_risk <- as.factor(mri_reduced$alcohol_risk)
mri_reduced$preference <- as.factor(mri_reduced$preference)
mri_reduced$with_food <- as.factor(mri_reduced$with_food)
mri_reduced$sex <- as.factor(mri_reduced$sex)
mri_reduced$ethnicity <- as.factor(mri_reduced$ethnicity)
mri_reduced$education <- as.factor(mri_reduced$education)
mri_reduced$smoking <- as.factor(mri_reduced$smoking)
mri_reduced$diet_quality <- as.factor(mri_reduced$diet_quality)
mri_reduced$sleep_quality <- as.factor(mri_reduced$sleep_quality)
mri_reduced$blood_pressure <- as.factor(mri_reduced$blood_pressure)
mri_reduced$BMI <- as.factor(mri_reduced$BMI)
mri_reduced$health_rating <- as.factor(mri_reduced$health_rating)
mri_reduced$ACDs <- ifelse(mri_reduced$ACDs == 2, 0, 1)
str(mri_reduced)




library(mediation)
# --- 设定变量 ---
exposure <- "alcohol_risk"
outcome <- "ACDs"
# 协变量列表 (基于你之前的反馈，自动兼容 _mer 后缀)
potential_covs <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                    "TDI", "smoking", "diet_quality", "sleep_quality",
                    "blood_pressure", "BMI", "health_rating", 
                    "triglycerides", "HDL", "LDL", "glucose", "urate", "p25781_i2")
# 中介变量列表
potential_meds <- c("p25013_i2","p25014_i2","p25020_i2", "p25016_i2",
                    "p25490_i2","p25491_i2","p25499_i2","p25507_i2","p25508_i2","p25513_i2",
                    "p25515_i2","p25516_i2","p25517_i2","p25518_i2",
                    "p25523_i2","p25524_i2","p25526_i2","p25527_i2","p25528_i2","p25529_i2",
                    "p25530_i2","p25532_i2","p25533_i2","p25531_i2","p25534_i2","p25535_i2",
                    "p25536_i2","p25537_i2","p25538_i2","p25539_i2","p25540_i2","p25541_i2")

# 2. 自动筛选数据框中确实存在的列
actual_covs <- intersect(potential_covs, names(mri))
actual_meds <- intersect(potential_meds, names(mri))
all_vars <- unique(c(exposure, "ACDs", actual_covs, actual_meds))
# 3. 剔除所有涉及变量的缺失值 (mediation分析对NA零容忍)
mri_final <- na.omit(mri_reduced[, all_vars])
# 4. 确保暴露变量是 Factor
mri_final[[exposure]] <- as.factor(mri_final[[exposure]])
cat("清洗完成。最终样本量：", nrow(mri_final), "\n")
cat("结局变量分布：\n"); print(table(mri_final$ACDs))

# 初始化存储结果的列表
results_list <- list()

for (med in actual_meds) {
  cat("\n正在处理中介变量:", med, "...")
  
  # 构建公式
  f_med <- as.formula(paste(med, "~", paste(c(exposure, actual_covs), collapse = " + ")))
  f_out <- as.formula(paste("ACDs ~", paste(c(exposure, med, actual_covs), collapse = " + ")))
  
  # 模型拟合
  fit_m <- lm(f_med, data = mri_final)
  fit_y <- glm(f_out, family = binomial(link = "logit"), data = mri_final)
  
  # 执行中介效应计算
  # 因为 alcohol_risk 有 3 个水平，必须指定 treat.value (如 "3") 和 control.value (如 "1")
  med_run <- tryCatch({
    mediate(fit_m, fit_y, 
            treat = exposure, 
            mediator = med,
            control.value = "1", # 请根据实际参考组修改
            treat.value = "3",   # 请根据实际高风险组修改
            boot = TRUE, 
            sims = 100)          # 正式发表建议改为 500 或 1000
  }, error = function(e) {
    cat(" 失败: ", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(med_run)) {
    # 提取关键指标
    res <- data.frame(
      Mediator = med,
      ACME_Est = med_run$d0,        # 间接效应
      ACME_P   = med_run$d0.p,
      ADE_Est  = med_run$z0,        # 直接效应
      ADE_P    = med_run$z0.p,
      Total_Est = med_run$tau.coef, # 总效应
      Total_P   = med_run$tau.p,
      Prop_Med  = med_run$n0,       # 中介占比
      stringsAsFactors = FALSE
    )
    results_list[[med]] <- res
    cat(" 成功!")
  }
}

# 合并所有结果
final_results <- do.call(rbind, results_list)
print(final_results)

# 保存结果
write.csv(final_results, "Mediation_Results_Alcohol_ACDs.csv", row.names = FALSE)




##### 补充分析
wmh_data <- read.csv("white matter hyperintensities.csv", stringsAsFactors = FALSE)
mri <- merge(mri, wmh_data[, c("eid", colnames(wmh_data)[2])], 
             by = "eid", all.x = TRUE)

library(boot)
library(dplyr)

# 1. 变量定义
fa_tracts <- c(
  "p25488_i2","p25489_i2","p25490_i2","p25491_i2",
  "p25492_i2","p25493_i2","p25496_i2","p25497_i2",
  "p25498_i2","p25499_i2","p25500_i2","p25501_i2",
  "p25502_i2","p25503_i2","p25505_i2","p25506_i2",
  "p25504_i2","p25494_i2","p25495_i2","p25507_i2",
  "p25508_i2","p25509_i2","p25510_i2","p25511_i2",
  "p25512_i2","p25513_i2","p25514_i2"
)

md_tracts <- c(
  "p25515_i2", "p25516_i2", "p25517_i2", "p25518_i2", 
  "p25519_i2", "p25520_i2", "p25523_i2", "p25524_i2", 
  "p25525_i2", "p25526_i2", "p25527_i2", "p25528_i2", 
  "p25529_i2", "p25530_i2", "p25532_i2", "p25533_i2", 
  "p25531_i2", "p25521_i2", "p25522_i2", "p25534_i2", 
  "p25535_i2", "p25536_i2", "p25537_i2", "p25538_i2", 
  "p25539_i2", "p25540_i2", "p25541_i2"
)

# 协变量向量（已删除末尾逗号）
covar_vec <- c(
  "preference", "with_food", "age", "sex", "ethnicity", 
  "education", "TDI", "smoking", "diet_quality", 
  "sleep_quality_mer", "blood_pressure_mer", "BMI", 
  "health_rating", "triglycerides", "HDL", "LDL", 
  "glucose", "urate"
)

library(boot)
library(dplyr)

# 1. 确保变量类型正确
mri$ACDs <- as.numeric(as.character(mri$ACDs))
mri$alcohol_risk <- as.factor(mri$alcohol_risk) # 确保是因子，且默认 1 为参考级

# 2. 定义计算函数
mediation_comparison_fn <- function(data, indices) {
  df_boot <- data[indices, ]
  
  calc_mean_mp <- function(tract_list, exposure_name, outcome, covars) {
    mps <- c()
    cov_str <- paste(covars, collapse = " + ")
    
    # 定义我们要提取的系数名称（以 1 为参考，提取 3 的系数）
    target_coef_exposure <- paste0(exposure_name, "3")
    
    for (m_name in tract_list) {
      mp_val <- tryCatch({
        # 模型 M: Mediator ~ Exposure + Covariates
        fm <- as.formula(paste(m_name, "~", exposure_name, "+", cov_str))
        fit_m <- lm(fm, data = df_boot)
        # 提取 alcohol_risk 等级 3 的系数 (a)
        a <- coef(fit_m)[target_coef_exposure]
        
        # 模型 Y: Outcome ~ Exposure + Mediator + Covariates (Probit)
        fy <- as.formula(paste(outcome, "~", exposure_name, "+", m_name, "+", cov_str))
        fit_y <- glm(fy, data = df_boot, family = binomial(link = "probit"))
        # 提取 Mediator 的系数 (b) 和 alcohol_risk 等级 3 的直接效应 (cp)
        b <- coef(fit_y)[m_name]
        cp <- coef(fit_y)[target_coef_exposure]
        
        indirect <- a * b
        total <- indirect + cp
        res <- indirect / total
        
        # 结果过滤：只保留有意义的比例（剔除极值）
        if(is.finite(res) && abs(res) <= 1) res else NA
      }, error = function(e) NA)
      
      if(!is.na(mp_val)) mps <- c(mps, mp_val)
    }
    
    if(length(mps) == 0) return(NA)
    return(mean(mps))
  }
  
  # 执行计算 (注意这里的 exposure_name 传入 "alcohol_risk")
  m_fa <- calc_mean_mp(fa_tracts, "alcohol_risk", "ACDs", covar_vec)
  m_md <- calc_mean_mp(md_tracts, "alcohol_risk", "ACDs", covar_vec)
  
  # 结果合并
  res_final <- c(
    if(is.na(m_fa)) 0 else m_fa,
    if(is.na(m_md)) 0 else m_md,
    if(is.na(m_fa) | is.na(m_md)) 0 else m_fa - m_md
  )
  return(res_final)
}

# 3. 运行自检
print("Running initial check with Factor adjustment...")
test_res <- mediation_comparison_fn(mri, 1:nrow(mri))
print("Test result (Mean FA MP, Mean MD MP, Difference):")
print(test_res)

# 4. 正式运行 Bootstrap
if(!all(test_res == 0)) {
  set.seed(123)
  boot_results <- boot(
    data = mri, 
    statistic = mediation_comparison_fn, 
    R = 1000, 
    parallel = "no" # 如果在服务器或非Windows系统，可改用 multicore
  )
  print(boot_results)
  
  # 提取 95% 置信区间
  ci_diff <- boot.ci(boot_results, type = "perc", index = 3)
  print(ci_diff)
}





##### 添加WMH作为协变量的中介分析
# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\alcohol投稿\\16-peerj major revision\\peerj major revision 新增材料")
getwd()
mri <- read.csv("mri_mediation with time.csv")
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
str(mri)
mri$alcohol_pw <- scale(mri$alcohol_pw)
## 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
## 由于由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并（将1和2合并为1-poor or moderate，3改为2-healthy）
## blood_pressure3、blood_pressure4的VIF大于5 需要进行成分合并（将1和2合并为1-normal or elevated，3改为2-stage 1, 4改为3-stage 2）
table(mri$blood_pressure)
chisq.test(table(mri$blood_pressure))
mri <- mri %>%
  mutate(blood_pressure_mer = case_when(
    blood_pressure %in% c(1, 2) ~ 1,
    blood_pressure == 3 ~ 2,
    blood_pressure == 4 ~ 3,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
mri$blood_pressure_mer <- as.factor(mri$blood_pressure_mer)
## 由于sleep_quality_2、sleep_quality_3的VIF大于5 需要进行成分合并
table(mri$sleep_quality)
chisq.test(table(mri$sleep_quality))
mri <- mri %>%
  mutate(sleep_quality_mer = case_when(
    sleep_quality %in% c(1, 2) ~ 1,
    sleep_quality == 3 ~ 2,
    TRUE ~ NA_real_ # 如果不符合上述任一条件，可以选择填入NA或适当的默认值
  )
  )
mri$sleep_quality_mer <- as.factor(mri$sleep_quality_mer)
wmh_data <- read.csv("white matter hyperintensities.csv", stringsAsFactors = FALSE)
mri <- merge(mri, wmh_data[, c("eid", colnames(wmh_data)[2])], 
             by = "eid", all.x = TRUE)
# 对 WMH 进行对数转换 (Field ID: p25781_i2)
mri$log_wmh <- log(mri$p25781_i2 + 1)
mri$ACDs[mri$ACDs == 2] <- 0
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


# 中介分析的敏感性分析
# 2. 更新协变量列表，确保名称与 mri 中的列名完全一致
exposure <- "alcohol_risk"
outcome  <- "ACDs"
# 注意：这里使用你最新合并生成的变量名
covariates <- c("preference", "with_food", "age", "sex", "ethnicity", "education", 
                "TDI", "smoking", "diet_quality", "sleep_quality_mer", "blood_pressure_mer", 
                "BMI", "health_rating", "triglycerides", "HDL", "LDL", "glucose", "urate", "log_wmh")

# 3. 确保中介变量列表正确
actual_meds <- c("p25013_i2","p25014_i2","p25020_i2","p25016_i2",
                 "p25490_i2","p25491_i2","p25499_i2","p25507_i2","p25508_i2","p25513_i2",
                 "p25515_i2","p25516_i2","p25517_i2","p25518_i2",
                 "p25523_i2","p25524_i2","p25526_i2","p25527_i2","p25528_i2","p25529_i2",
                 "p25530_i2","p25532_i2","p25533_i2","p25531_i2","p25534_i2","p25535_i2",
                 "p25536_i2","p25537_i2","p25538_i2","p25539_i2","p25540_i2","p25541_i2")

# 4. 执行循环
results_list <- list()

for (med in actual_meds) {
  cat("\n正在处理中介变量:", med, "...")
  
  # 检查中介变量是否存在于 mri 中
  if (!(med %in% names(mri))) {
    cat(" 变量不存在于 mri 中，跳过。")
    next
  }
  
  # 构建公式
  f_med <- as.formula(paste(med, "~", exposure, "+", paste(covariates, collapse = " + ")))
  f_out <- as.formula(paste(outcome, "~", exposure, "+", med, "+", paste(covariates, collapse = " + ")))
  
  med_run <- tryCatch({
    # 使用 mri 数据框进行拟合
    fit_m <- lm(f_med, data = mri)
    fit_y <- glm(f_out, family = binomial(link = "logit"), data = mri)
    
    # 执行中介分析 (boot=FALSE 使用准则法，大样本更稳健)
    mediate(fit_m, fit_y, 
            treat = exposure, 
            mediator = med,
            control.value = "1", 
            treat.value = "3",
            boot = FALSE,   
            sims = 500)     
  }, error = function(e) {
    cat(" 失败: ", e$message)
    return(NULL)
  })
  
  if (!is.null(med_run)) {
    # 提取结果
    res <- data.frame(
      Mediator = med,
      ACME_Est = med_run$d0, 
      ACME_P   = med_run$d0.p,
      ADE_Est  = med_run$z0, 
      ADE_P    = med_run$z0.p,
      Total_Est = med_run$tau.coef,
      Prop_Med = med_run$n0,
      stringsAsFactors = FALSE
    )
    results_list[[med]] <- res
    cat(" 成功!")
  }
}

# 汇总结果查看
final_results <- do.call(rbind, results_list)
print(final_results)