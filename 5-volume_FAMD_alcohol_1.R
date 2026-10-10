# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mri_mediation.csv")

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
mri$ACDs <- as.factor(mri$ACDs)
mri$All_cause_stroke <- as.factor(mri$All_cause_stroke)
mri$Hemorrhagic_stroke <- as.factor(mri$Hemorrhagic_stroke)
mri$Ischemic_stroke <- as.factor(mri$Ischemic_stroke)
mri$TIA <- as.factor(mri$TIA)
str(mri)
mri$alcohol_pw <- scale(mri$alcohol_pw)


##### 模型1
##### 皮层下区域

# 检查模型1协变量之间的多重共线性
model_1 <- vif(lm(alcohol_pw ~ preference + with_food, data = mri))
print(model_1)
# accumbens_left
model1_1L <- lm(p25023_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_1L)
# accumbens_right
model1_1R <- lm(p25024_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_1R)
# amygdala_left
model1_2L <- lm(p25021_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_2L)
# amygdala_right
model1_2R <- lm(p25022_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_2R)
# caudate_left
model1_3L <- lm(p25013_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_3L)
# caudate_right
model1_3R <- lm(p25014_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_3R)
# hippocampus_left
model1_4L <- lm(p25019_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_4L)
# hippocampus_right
model1_4R <- lm(p25020_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_4R)
# pallidum_left
model1_5L <- lm(p25017_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_5L)
# pallidum_right
model1_5R <- lm(p25018_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_5R)
# putamen_left
model1_6L <- lm(p25015_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_6L)
# putamen_right
model1_6R <- lm(p25016_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_6R)
# thalamus_left
model1_7L <- lm(p25011_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_7L)
# thalamus_right
model1_7R <- lm(p25012_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_7R)


##### 各个脑区的体积
# cuneus_left
model1_8L <- lm(p27207_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_8L)
# cuneus_right
model1_8R <- lm(p27300_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_8R)
# entorhinal_left
model1_9L <- lm(p27208_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_9L)
# entorhinal_right
model1_9R <- lm(p27301_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_9R)
# fusiform_left
model1_10L <- lm(p27209_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_10L)
# fusiform_right
model1_10R <- lm(p27302_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_10R)
# inferiorparietal_left
model1_11L <- lm(p27210_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_11L)
# inferiorparietal_right
model1_11R <- lm(p27303_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_11R)
# inferiortemporal_left
model1_12L <- lm(p27211_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_12L)
# inferiortemporal_right
model1_12R <- lm(p27304_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_12R)
# insula_left
model1_13L <- lm(p27235_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_13L)
# insula_right
model1_13R <- lm(p27328_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_13R)
# isthmuscingulate_left
model1_14L <- lm(p27212_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_14L)
# isthmuscingulate_right
model1_14R <- lm(p27305_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_14R)
# lateraloccipital_left
model1_15L <- lm(p27213_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_15L)
# lateraloccipital_right
model1_15R <- lm(p27306_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_15R)
# lateralorbitofrontal_left
model1_16L <- lm(p27214_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_16L)
# lateralorbitofrontal_right
model1_16R <- lm(p27307_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_16R)
# lingual_left
model1_17L <- lm(p27215_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_17L)
# lingual_right
model1_17R <- lm(p27308_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_17R)
# medialorbitofrontal_left
model1_18L <- lm(p27216_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_18L)
# medialorbitofrontal_right
model1_18R <- lm(p27309_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_18R)
# middletemporal_left
model1_19L <- lm(p27217_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_19L)
# middletemporal_right
model1_19R <- lm(p27310_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_19R)
# paracentral_left
model1_20L <- lm(p27219_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_20L)
# paracentral_right
model1_20R <- lm(p27312_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_20R)
# parahippocampal_left
model1_21L <- lm(p27218_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_21L)
# parahippocampal_right
model1_21R <- lm(p27311_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_21R)
# parsopercularis_left
model1_22L <- lm(p27220_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_22L)
# parsopercularis_right
model1_22R <- lm(p27313_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_22R)
# parsorbitalis_left
model1_23L <- lm(p27221_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_23L)
# parsorbitalis_right
model1_23R <- lm(p27314_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_23R)
# parstriangularis_left
model1_24L <- lm(p27222_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_24L)
# parstriangularis_right
model1_24R <- lm(p27315_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_24R)
# pericalcarine_left
model1_25L <- lm(p27223_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_25L)
# pericalcarine_right
model1_25R <- lm(p27316_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_25R)
# postcentral_left
model1_26L <- lm(p27224_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_26L)
# postcentral_right
model1_26R <- lm(p27317_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_26R)
# posteriorcingulate_left
model1_27L <- lm(p27225_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_27L)
# posteriorcingulate_right
model1_27R <- lm(p27318_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_27R)
# precentral_left
model1_28L <- lm(p27226_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_28L)
# precentral_right
model1_28R <- lm(p27319_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_28R)
# precuneus_left
model1_29L <- lm(p27227_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_29L)
# precuneus_right
model1_29R <- lm(p27320_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_29R)
# rostralanteriorcingulate_left
model1_30L <- lm(p27228_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_30L)
# rostralanteriorcingulate_right
model1_30R <- lm(p27321_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_30R)
# rostralmiddlefrontal_left
model1_31L <- lm(p27229_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_31L)
# rostralmiddlefrontal_right
model1_31R <- lm(p27322_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_31R)
# superiorfrontal_left
model1_32L <- lm(p27230_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_32L)
# superiorfrontal_right
model1_32R <- lm(p27323_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_32R)
# superiorparietal_left
model1_33L<- lm(p27231_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_33L)
# superiorparietal_right
model1_33R <- lm(p27324_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_33R)
# superiortemporal_left
model1_34L <- lm(p27232_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_34L)
# superiortemporal_right
model1_34R <- lm(p27325_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_34R)
# supramarginal_left
model1_35L <- lm(p27233_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_35L)
# supramarginal_right
model1_35R <- lm(p27326_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_35R)
# transversetemporal_left
model1_36L <- lm(p27234_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_36L)
# transversetemporal_right
model1_36R <- lm(p27327_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_36R)
# caudalanteriorcingulate_left
model1_37L <- lm(p27205_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_37L)
# caudalanteriorcingulate_right
model1_37R <- lm(p27298_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_37R)
# caudalmiddlefrontal_left
model1_38L <- lm(p27206_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_38L)
# caudalmiddlefrontal_right
model1_38R <- lm(p27299_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_38R)
# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model1_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_left
  summary(model1_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_right
  summary(model1_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_left
  summary(model1_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_right
  summary(model1_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_left
  summary(model1_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_right
  summary(model1_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_left
  summary(model1_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_right
  summary(model1_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_left
  summary(model1_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_right
  summary(model1_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_left
  summary(model1_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_right
  summary(model1_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_left
  summary(model1_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_right
  summary(model1_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_left
  summary(model1_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_right
  summary(model1_9L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_left
  summary(model1_9R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_right
  summary(model1_10L)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_left
  summary(model1_10R)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_right
  summary(model1_11L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_left
  summary(model1_11R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_right
  summary(model1_12L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_left
  summary(model1_12R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_right
  summary(model1_13L)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_left
  summary(model1_13R)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_right
  summary(model1_14L)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_left
  summary(model1_14R)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_right
  summary(model1_15L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_left
  summary(model1_15R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_right
  summary(model1_16L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_left
  summary(model1_16R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_right
  summary(model1_17L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_left
  summary(model1_17R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_right
  summary(model1_18L)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_left
  summary(model1_18R)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_right
  summary(model1_19L)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_left
  summary(model1_19R)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_right
  summary(model1_20L)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_left
  summary(model1_20R)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_right
  summary(model1_21L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_left
  summary(model1_21R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_right
  summary(model1_22L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_left
  summary(model1_22R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_right
  summary(model1_23L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_left
  summary(model1_23R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_right
  summary(model1_24L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_left
  summary(model1_24R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_right
  summary(model1_25L)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_left
  summary(model1_25R)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_right
  summary(model1_26L)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_left
  summary(model1_26R)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_right
  summary(model1_27L)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_left
  summary(model1_27R)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_right
  summary(model1_28L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_left
  summary(model1_28R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_right
  summary(model1_29L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_left
  summary(model1_29R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_right
  summary(model1_30L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_left
  summary(model1_30R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_right
  summary(model1_31L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_left
  summary(model1_31R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_right
  summary(model1_32L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_left
  summary(model1_32R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_right
  summary(model1_33L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_left
  summary(model1_33R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_right
  summary(model1_34L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_left
  summary(model1_34R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_right
  summary(model1_35L)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_left
  summary(model1_35R)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_right
  summary(model1_36L)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_left
  summary(model1_36R)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_right
  summary(model1_37L)$coefficients["alcohol_pw", "Pr(>|t|)"], # caudalanteriorcingulate_left
  summary(model1_37R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalanteriorcingulate_right
  summary(model1_38L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalmiddlefrontal_left
  summary(model1_38R)$coefficients["alcohol_pw", "Pr(>|t|)"]   # caudalmiddlefrontal_right
)

# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)

##### 模型2
##### 皮层下区域

# 检查模型2协变量的多重线性关系
model_2 <- vif(lm(alcohol_pw ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri))
print(model_2)
# accumbens_left
model2_1L <- lm(p25023_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_1L)
# accumbens_right
model2_1R <- lm(p25024_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_1R)
# amygdala_left
model2_2L <- lm(p25021_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_2L)
# amygdala_right
model2_2R <- lm(p25022_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_2R)
# caudate_left
model2_3L <- lm(p25013_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_3L)
# caudate_right
model2_3R <- lm(p25014_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_3R)
# hippocampus_left
model2_4L <- lm(p25019_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_4L)
# hippocampus_right
model2_4R <- lm(p25020_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_4R)
# pallidum_left
model2_5L <- lm(p25017_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_5L)
# pallidum_right
model2_5R <- lm(p25018_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_5R)
# putamen_left
model2_6L <- lm(p25015_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_6L)
# putamen_right
model2_6R <- lm(p25016_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_6R)
# thalamus_left
model2_7L <- lm(p25011_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_7L)
# thalamus_right
model2_7R <- lm(p25012_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_7R)


##### 各个脑区的体积
# cuneus_left
model2_8L <- lm(p27207_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_8L)
# cuneus_right
model2_8R <- lm(p27300_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_8R)
# entorhinal_left
model2_9L <- lm(p27208_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_9L)
# entorhinal_right
model2_9R <- lm(p27301_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_9R)
# fusiform_left
model2_10L <- lm(p27209_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_10L)
# fusiform_right
model2_10R <- lm(p27302_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_10R)
# inferiorparietal_left
model2_11L <- lm(p27210_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_11L)
# inferiorparietal_right
model2_11R <- lm(p27303_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_11R)
# inferiortemporal_left
model2_12L <- lm(p27211_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_12L)
# inferiortemporal_right
model2_12R <- lm(p27304_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_12R)
# insula_left
model2_13L <- lm(p27235_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_13L)
# insula_right
model2_13R <- lm(p27328_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_13R)
# isthmuscingulate_left
model2_14L <- lm(p27212_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_14L)
# isthmuscingulate_right
model2_14R <- lm(p27305_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_14R)
# lateraloccipital_left
model2_15L <- lm(p27213_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_15L)
# lateraloccipital_right
model2_15R <- lm(p27306_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_15R)
# lateralorbitofrontal_left
model2_16L <- lm(p27214_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_16L)
# lateralorbitofrontal_right
model2_16R <- lm(p27307_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_16R)
# lingual_left
model2_17L <- lm(p27215_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_17L)
# lingual_right
model2_17R <- lm(p27308_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_17R)
# medialorbitofrontal_left
model2_18L <- lm(p27216_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_18L)
# medialorbitofrontal_right
model2_18R <- lm(p27309_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_18R)
# middletemporal_left
model2_19L <- lm(p27217_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_19L)
# middletemporal_right
model2_19R <- lm(p27310_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_19R)
# paracentral_left
model2_20L <- lm(p27219_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_20L)
# paracentral_right
model2_20R <- lm(p27312_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_20R)
# parahippocampal_left
model2_21L <- lm(p27218_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_21L)
# parahippocampal_right
model2_21R <- lm(p27311_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_21R)
# parsopercularis_left
model2_22L <- lm(p27220_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_22L)
# parsopercularis_right
model2_22R <- lm(p27313_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_22R)
# parsorbitalis_left
model2_23L <- lm(p27221_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_23L)
# parsorbitalis_right
model2_23R <- lm(p27314_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_23R)
# parstriangularis_left
model2_24L <- lm(p27222_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_24L)
# parstriangularis_right
model2_24R <- lm(p27315_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_24R)
# pericalcarine_left
model2_25L <- lm(p27223_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_25L)
# pericalcarine_right
model2_25R <- lm(p27316_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_25R)
# postcentral_left
model2_26L <- lm(p27224_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_26L)
# postcentral_right
model2_26R <- lm(p27317_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_26R)
# posteriorcingulate_left
model2_27L <- lm(p27225_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_27L)
# posteriorcingulate_right
model2_27R <- lm(p27318_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_27R)
# precentral_left
model2_28L <- lm(p27226_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_28L)
# precentral_right
model2_28R <- lm(p27319_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_28R)
# precuneus_left
model2_29L <- lm(p27227_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_29L)
# precuneus_right
model2_29R <- lm(p27320_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_29R)
# rostralanteriorcingulate_left
model2_30L <- lm(p27228_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_30L)
# rostralanteriorcingulate_right
model2_30R <- lm(p27321_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_30R)
# rostralmiddlefrontal_left
model2_31L <- lm(p27229_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_31L)
# rostralmiddlefrontal_right
model2_31R <- lm(p27322_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_31R)
# superiorfrontal_left
model2_32L <- lm(p27230_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_32L)
# superiorfrontal_right
model2_32R <- lm(p27323_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_32R)
# superiorparietal_left
model2_33L<- lm(p27231_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_33L)
# superiorparietal_right
model2_33R <- lm(p27324_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_33R)
# superiortemporal_left
model2_34L <- lm(p27232_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_34L)
# superiortemporal_right
model2_34R <- lm(p27325_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_34R)
# supramarginal_left
model2_35L <- lm(p27233_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_35L)
# supramarginal_right
model2_35R <- lm(p27326_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_35R)
# transversetemporal_left
model2_36L <- lm(p27234_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_36L)
# transversetemporal_right
model2_36R <- lm(p27327_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_36R)
# caudalanteriorcingulate_left
model2_37L <- lm(p27205_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_37L)
# caudalanteriorcingulate_right
model2_37R <- lm(p27298_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_37R)
# caudalmiddlefrontal_left
model2_38L <- lm(p27206_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_38L)
# caudalmiddlefrontal_right
model2_38R <- lm(p27299_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_38R)
# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model2_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_left
  summary(model2_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_right
  summary(model2_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_left
  summary(model2_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_right
  summary(model2_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_left
  summary(model2_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_right
  summary(model2_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_left
  summary(model2_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_right
  summary(model2_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_left
  summary(model2_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_right
  summary(model2_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_left
  summary(model2_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_right
  summary(model2_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_left
  summary(model2_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_right
  summary(model2_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_left
  summary(model2_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_right
  summary(model2_9L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_left
  summary(model2_9R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_right
  summary(model2_10L)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_left
  summary(model2_10R)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_right
  summary(model2_11L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_left
  summary(model2_11R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_right
  summary(model2_12L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_left
  summary(model2_12R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_right
  summary(model2_13L)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_left
  summary(model2_13R)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_right
  summary(model2_14L)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_left
  summary(model2_14R)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_right
  summary(model2_15L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_left
  summary(model2_15R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_right
  summary(model2_16L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_left
  summary(model2_16R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_right
  summary(model2_17L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_left
  summary(model2_17R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_right
  summary(model2_18L)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_left
  summary(model2_18R)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_right
  summary(model2_19L)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_left
  summary(model2_19R)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_right
  summary(model2_20L)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_left
  summary(model2_20R)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_right
  summary(model2_21L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_left
  summary(model2_21R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_right
  summary(model2_22L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_left
  summary(model2_22R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_right
  summary(model2_23L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_left
  summary(model2_23R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_right
  summary(model2_24L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_left
  summary(model2_24R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_right
  summary(model2_25L)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_left
  summary(model2_25R)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_right
  summary(model2_26L)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_left
  summary(model2_26R)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_right
  summary(model2_27L)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_left
  summary(model2_27R)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_right
  summary(model2_28L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_left
  summary(model2_28R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_right
  summary(model2_29L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_left
  summary(model2_29R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_right
  summary(model2_30L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_left
  summary(model2_30R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_right
  summary(model2_31L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_left
  summary(model2_31R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_right
  summary(model2_32L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_left
  summary(model2_32R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_right
  summary(model2_33L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_left
  summary(model2_33R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_right
  summary(model2_34L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_left
  summary(model2_34R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_right
  summary(model2_35L)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_left
  summary(model2_35R)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_right
  summary(model2_36L)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_left
  summary(model2_36R)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_right
  summary(model2_37L)$coefficients["alcohol_pw", "Pr(>|t|)"], # caudalanteriorcingulate_left
  summary(model2_37R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalanteriorcingulate_right
  summary(model2_38L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalmiddlefrontal_left
  summary(model2_38R)$coefficients["alcohol_pw", "Pr(>|t|)"]   # caudalmiddlefrontal_right
)

# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)

##### 模型3
##### 皮层下区域

# 检查模型3协变量之间的多重共线性
model_3 <- vif(lm(alcohol_pw ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + cholesterol + triglycerides + HDL + LDL + glucose + urate , data = mri))
print(model_3)
# 由于cholesterol、triglycerides、HDL、LDL的VIF较高，因此去掉cholesterol进行分析
model_3 <- vif(lm(alcohol_pw ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = mri))
print(model_3)

# accumbens_left
model3_1L <- lm(p25023_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_1L)
# accumbens_right
model3_1R <- lm(p25024_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_1R)
# amygdala_left
model3_2L <- lm(p25021_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_2L)
# amygdala_right
model3_2R <- lm(p25022_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_2R)
# caudate_left
model3_3L <- lm(p25013_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_3L)
# caudate_right
model3_3R <- lm(p25014_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_3R)
# hippocampus_left
model3_4L <- lm(p25019_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_4L)
# hippocampus_right
model3_4R <- lm(p25020_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_4R)
# pallidum_left
model3_5L <- lm(p25017_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_5L)
# pallidum_right
model3_5R <- lm(p25018_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_5R)
# putamen_left
model3_6L <- lm(p25015_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_6L)
# putamen_right
model3_6R <- lm(p25016_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_6R)
# thalamus_left
model3_7L <- lm(p25011_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_7L)
# thalamus_right
model3_7R <- lm(p25012_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_7R)

##### 各个脑区的体积
# cuneus_left
model3_8L <- lm(p27207_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_8L)
# cuneus_right
model3_8R <- lm(p27300_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_8R)
# entorhinal_left
model3_9L <- lm(p27208_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_9L)
# entorhinal_right
model3_9R <- lm(p27301_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_9R)
# fusiform_left
model3_10L <- lm(p27209_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_10L)
# fusiform_right
model3_10R <- lm(p27302_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_10R)
# inferiorparietal_left
model3_11L <- lm(p27210_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_11L)
# inferiorparietal_right
model3_11R <- lm(p27303_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_11R)
# inferiortemporal_left
model3_12L <- lm(p27211_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_12L)
# inferiortemporal_right
model3_12R <- lm(p27304_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_12R)
# insula_left
model3_13L <- lm(p27235_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_13L)
# insula_right
model3_13R <- lm(p27328_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_13R)
# isthmuscingulate_left
model3_14L <- lm(p27212_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_14L)
# isthmuscingulate_right
model3_14R <- lm(p27305_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_14R)
# lateraloccipital_left
model3_15L <- lm(p27213_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_15L)
# lateraloccipital_right
model3_15R <- lm(p27306_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_15R)
# lateralorbitofrontal_left
model3_16L <- lm(p27214_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_16L)
# lateralorbitofrontal_right
model3_16R <- lm(p27307_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_16R)
# lingual_left
model3_17L <- lm(p27215_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_17L)
# lingual_right
model3_17R <- lm(p27308_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_17R)
# medialorbitofrontal_left
model3_18L <- lm(p27216_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_18L)
# medialorbitofrontal_right
model3_18R <- lm(p27309_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_18R)
# middletemporal_left
model3_19L <- lm(p27217_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_19L)
# middletemporal_right
model3_19R <- lm(p27310_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_19R)
# paracentral_left
model3_20L <- lm(p27219_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_20L)
# paracentral_right
model3_20R <- lm(p27312_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_20R)
# parahippocampal_left
model3_21L <- lm(p27218_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_21L)
# parahippocampal_right
model3_21R <- lm(p27311_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_21R)
# parsopercularis_left
model3_22L <- lm(p27220_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_22L)
# parsopercularis_right
model3_22R <- lm(p27313_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_22R)
# parsorbitalis_left
model3_23L <- lm(p27221_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_23L)
# parsorbitalis_right
model3_23R <- lm(p27314_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_23R)
# parstriangularis_left
model3_24L <- lm(p27222_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_24L)
# parstriangularis_right
model3_24R <- lm(p27315_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_24R)
# pericalcarine_left
model3_25L <- lm(p27223_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_25L)
# pericalcarine_right
model3_25R <- lm(p27316_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_25R)
# postcentral_left
model3_26L <- lm(p27224_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_26L)
# postcentral_right
model3_26R <- lm(p27317_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_26R)
# posteriorcingulate_left
model3_27L <- lm(p27225_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_27L)
# posteriorcingulate_right
model3_27R <- lm(p27318_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_27R)
# precentral_left
model3_28L <- lm(p27226_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_28L)
# precentral_right
model3_28R <- lm(p27319_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_28R)
# precuneus_left
model3_29L <- lm(p27227_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_29L)
# precuneus_right
model3_29R <- lm(p27320_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_29R)
# rostralanteriorcingulate_left
model3_30L <- lm(p27228_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_30L)
# rostralanteriorcingulate_right
model3_30R <- lm(p27321_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_30R)
# rostralmiddlefrontal_left
model3_31L <- lm(p27229_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating+ triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_31L)
# rostralmiddlefrontal_right
model3_31R <- lm(p27322_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_31R)
# superiorfrontal_left
model3_32L <- lm(p27230_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating+ triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_32L)
# superiorfrontal_right
model3_32R <- lm(p27323_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating+ triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_32R)
# superiorparietal_left
model3_33L<- lm(p27231_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_33L)
# superiorparietal_right
model3_33R <- lm(p27324_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_33R)
# superiortemporal_left
model3_34L <- lm(p27232_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_34L)
# superiortemporal_right
model3_34R <- lm(p27325_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_34R)
# supramarginal_left
model3_35L <- lm(p27233_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_35L)
# supramarginal_right
model3_35R <- lm(p27326_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_35R)
# transversetemporal_left
model3_36L <- lm(p27234_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_36L)
# transversetemporal_right
model3_36R <- lm(p27327_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_36R)
# caudalanteriorcingulate_left
model3_37L <- lm(p27205_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_37L)
# caudalanteriorcingulate_right
model3_37R <- lm(p27298_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_37R)
# caudalmiddlefrontal_left
model3_38L <- lm(p27206_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_38L)
# caudalmiddlefrontal_right
model3_38R <- lm(p27299_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_38R)

# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model3_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_left
  summary(model3_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # accumbens_right
  summary(model3_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_left
  summary(model3_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # amygdala_right
  summary(model3_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_left
  summary(model3_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudate_right
  summary(model3_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_left
  summary(model3_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # hippocampus_right
  summary(model3_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_left
  summary(model3_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # pallidum_right
  summary(model3_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_left
  summary(model3_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # putamen_right
  summary(model3_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_left
  summary(model3_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # thalamus_right
  summary(model3_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_left
  summary(model3_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # cuneus_right
  summary(model3_9L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_left
  summary(model3_9R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # entorhinal_right
  summary(model3_10L)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_left
  summary(model3_10R)$coefficients["alcohol_pw", "Pr(>|t|)"], # fusiform_right
  summary(model3_11L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_left
  summary(model3_11R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiorparietal_right
  summary(model3_12L)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_left
  summary(model3_12R)$coefficients["alcohol_pw", "Pr(>|t|)"], # inferiortemporal_right
  summary(model3_13L)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_left
  summary(model3_13R)$coefficients["alcohol_pw", "Pr(>|t|)"], # insula_right
  summary(model3_14L)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_left
  summary(model3_14R)$coefficients["alcohol_pw", "Pr(>|t|)"], # isthmuscingulate_right
  summary(model3_15L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_left
  summary(model3_15R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateraloccipital_right
  summary(model3_16L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_left
  summary(model3_16R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lateralorbitofrontal_right
  summary(model3_17L)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_left
  summary(model3_17R)$coefficients["alcohol_pw", "Pr(>|t|)"], # lingual_right
  summary(model3_18L)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_left
  summary(model3_18R)$coefficients["alcohol_pw", "Pr(>|t|)"], # medialorbitofrontal_right
  summary(model3_19L)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_left
  summary(model3_19R)$coefficients["alcohol_pw", "Pr(>|t|)"], # middletemporal_right
  summary(model3_20L)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_left
  summary(model3_20R)$coefficients["alcohol_pw", "Pr(>|t|)"], # paracentral_right
  summary(model3_21L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_left
  summary(model3_21R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parahippocampal_right
  summary(model3_22L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_left
  summary(model3_22R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsopercularis_right
  summary(model3_23L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_left
  summary(model3_23R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parsorbitalis_right
  summary(model3_24L)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_left
  summary(model3_24R)$coefficients["alcohol_pw", "Pr(>|t|)"], # parstriangularis_right
  summary(model3_25L)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_left
  summary(model3_25R)$coefficients["alcohol_pw", "Pr(>|t|)"], # pericalcarine_right
  summary(model3_26L)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_left
  summary(model3_26R)$coefficients["alcohol_pw", "Pr(>|t|)"], # postcentral_right
  summary(model3_27L)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_left
  summary(model3_27R)$coefficients["alcohol_pw", "Pr(>|t|)"], # posteriorcingulate_right
  summary(model3_28L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_left
  summary(model3_28R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precentral_right
  summary(model3_29L)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_left
  summary(model3_29R)$coefficients["alcohol_pw", "Pr(>|t|)"], # precuneus_right
  summary(model3_30L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_left
  summary(model3_30R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralanteriorcingulate_right
  summary(model3_31L)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_left
  summary(model3_31R)$coefficients["alcohol_pw", "Pr(>|t|)"], # rostralmiddlefrontal_right
  summary(model3_32L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_left
  summary(model3_32R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorfrontal_right
  summary(model3_33L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_left
  summary(model3_33R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiorparietal_right
  summary(model3_34L)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_left
  summary(model3_34R)$coefficients["alcohol_pw", "Pr(>|t|)"], # superiortemporal_right
  summary(model3_35L)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_left
  summary(model3_35R)$coefficients["alcohol_pw", "Pr(>|t|)"], # supramarginal_right
  summary(model3_36L)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_left
  summary(model3_36R)$coefficients["alcohol_pw", "Pr(>|t|)"], # transversetemporal_right
  summary(model3_37L)$coefficients["alcohol_pw", "Pr(>|t|)"], # caudalanteriorcingulate_left
  summary(model3_37R)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalanteriorcingulate_right
  summary(model3_38L)$coefficients["alcohol_pw", "Pr(>|t|)"],  # caudalmiddlefrontal_left
  summary(model3_38R)$coefficients["alcohol_pw", "Pr(>|t|)"]   # caudalmiddlefrontal_right
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
model1_1L <- lm(p25488_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_1L)
# acoustic radiation_right
model1_1R <- lm(p25489_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_1R)
# anterior thalamic radiation_left
model1_2L <- lm(p25490_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_2L)
# anterior thalamic radiation_right
model1_2R <- lm(p25491_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_2R)
# cingulate gyrus part of cingulum_left
model1_3L <- lm(p25492_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_3L)
# cingulate gyrus part of cingulum_right
model1_3R <- lm(p25493_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_3R)
# corticospinal tract_left
model1_4L <- lm(p25496_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_4L)
# corticospinal tract_right
model1_4R <- lm(p25497_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_4R)
# forceps major_left
model1_5L <- lm(p25498_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_5L)
# forceps major_right
model1_5R <- lm(p25499_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_5R)
# inferior fronto-occipital fasciculus_left
model1_6L <- lm(p25500_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_6L)
# inferior fronto-occipital fasciculus_right
model1_6R <- lm(p25501_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_6R)
# longitudinal fasciculus_left
model1_7L <- lm(p25502_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_7L)
# longitudinal fasciculus_right
model1_7R <- lm(p25503_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_7R)
# medial lemniscus_left
model1_8L <- lm(p25505_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_8L)
# medial lemniscus_right
model1_8R <- lm(p25506_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_8R)
# middle cerebellar peduncle
model1_9 <- lm(p25504_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_9)
# parahippocampal part of cingulum_left
model1_10L <- lm(p25494_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_10L)
# parahippocampal part of cingulum_right
model1_10R <- lm(p25495_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_10R)
# posterior thalamic radiation_left
model1_11L <- lm(p25507_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_11L)
# posterior thalamic radiation_right
model1_11R <- lm(p25508_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_11R)
# superior longitudinal fasciculus_left
model1_12L <- lm(p25509_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_12L)
# superior longitudinal fasciculus_right
model1_12R <- lm(p25510_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_12R)
# superior thalamic radiation_left
model1_13L <- lm(p25511_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_13L)
# superior thalamic radiation_right
model1_13R <- lm(p25512_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_13R)
# uncinate fasciculus_left
model1_14L <- lm(p25513_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_14L)
# uncinate fasciculus_right
model1_14R <- lm(p25514_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_14R)
##### 模型1
##### MD
# acoustic radiation_left
model1_15L <- lm(p25515_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_15L)
# acoustic radiation_right
model1_15R <- lm(p25516_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_15R)
# anterior thalamic radiation_left
model1_16L <- lm(p25517_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_16L)
# anterior thalamic radiation_right
model1_16R <- lm(p25518_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_16R)
# cingulate gyrus part of cingulum_left
model1_17L <- lm(p25519_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_17L)
# cingulate gyrus part of cingulum_right
model1_17R <- lm(p25520_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_17R)
# corticospinal tract_left
model1_18L <- lm(p25523_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_18L)
# corticospinal tract_right
model1_18R <- lm(p25524_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_18R)
# forceps major_left
model1_19L <- lm(p25525_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_19L)
# forceps major_right
model1_19R <- lm(p25526_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_19R)
# inferior fronto-occipital fasciculus_left
model1_20L <- lm(p25527_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_20L)
# inferior fronto-occipital fasciculus_right
model1_20R <- lm(p25528_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_20R)
# longitudinal fasciculus_left
model1_21L <- lm(p25529_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_21L)
# longitudinal fasciculus_right
model1_21R <- lm(p25530_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_21R)
# medial lemniscus_left
model1_22L <- lm(p25532_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_22L)
# medial lemniscus_right
model1_22R <- lm(p25533_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_22R)
# middle cerebellar peduncle
model1_23 <- lm(p25531_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_23)
# parahippocampal part of cingulum_left
model1_24L <- lm(p25521_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_24L)
# parahippocampal part of cingulum_right
model1_24R <- lm(p25522_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_24R)
# posterior thalamic radiation_left
model1_25L <- lm(p25534_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_25L)
# posterior thalamic radiation_right
model1_25R <- lm(p25535_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_25R)
# superior longitudinal fasciculus_left
model1_26L <- lm(p25536_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_26L)
# superior longitudinal fasciculus_right
model1_26R <- lm(p25537_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_26R)
# superior thalamic radiation_left
model1_27L <- lm(p25538_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_27L)
# superior thalamic radiation_right
model1_27R <- lm(p25539_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_27R)
# uncinate fasciculus_left
model1_28L <- lm(p25540_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_28L)
# uncinate fasciculus_right
model1_28R <- lm(p25541_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model1_28R)
# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model1_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  
  summary(model1_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_9)$coefficients["alcohol_pw", "Pr(>|t|)"], 
  summary(model1_10L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_10R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_11L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_11R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_12L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_12R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_13L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_13R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_14L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_14R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_15L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_15R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_16L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_16R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_17L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_17R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_18L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_18R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_19L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_19R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_20L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_20R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_21L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_21R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_22L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_22R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_23)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_24L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_24R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_25L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_25R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_26L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_26R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_27L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_27R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_28L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model1_28R)$coefficients["alcohol_pw", "Pr(>|t|)"]
)

# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)

##### 模型2
##### FA
# acoustic radiation_left
model2_1L <- lm(p25488_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_1L)
# acoustic radiation_right
model2_1R <- lm(p25489_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_1R)
# anterior thalamic radiation_left
model2_2L <- lm(p25490_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_2L)
# anterior thalamic radiation_right
model2_2R <- lm(p25491_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_2R)
# cingulate gyrus part of cingulum_left
model2_3L <- lm(p25492_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_3L)
# cingulate gyrus part of cingulum_right
model2_3R <- lm(p25493_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_3R)
# corticospinal tract_left
model2_4L <- lm(p25496_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_4L)
# corticospinal tract_right
model2_4R <- lm(p25497_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_4R)
# forceps major_left
model2_5L <- lm(p25498_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_5L)
# forceps major_right
model2_5R <- lm(p25499_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_5R)
# inferior fronto-occipital fasciculus_left
model2_6L <- lm(p25500_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_6L)
# inferior fronto-occipital fasciculus_right
model2_6R <- lm(p25501_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_6R)
# longitudinal fasciculus_left
model2_7L <- lm(p25502_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_7L)
# longitudinal fasciculus_right
model2_7R <- lm(p25503_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_7R)
# medial lemniscus_left
model2_8L <- lm(p25505_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_8L)
# medial lemniscus_right
model2_8R <- lm(p25506_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_8R)
# middle cerebellar peduncle
model2_9 <- lm(p25504_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_9)
# parahippocampal part of cingulum_left
model2_10L <- lm(p25494_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_10L)
# parahippocampal part of cingulum_right
model2_10R <- lm(p25495_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_10R)
# posterior thalamic radiation_left
model2_11L <- lm(p25507_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_11L)
# posterior thalamic radiation_right
model2_11R <- lm(p25508_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_11R)
# superior longitudinal fasciculus_left
model2_12L <- lm(p25509_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_12L)
# superior longitudinal fasciculus_right
model2_12R <- lm(p25510_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_12R)
# superior thalamic radiation_left
model2_13L <- lm(p25511_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_13L)
# superior thalamic radiation_right
model2_13R <- lm(p25512_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_13R)
# uncinate fasciculus_left
model2_14L <- lm(p25513_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_14L)
# uncinate fasciculus_right
model2_14R <- lm(p25514_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_14R)
##### 模型2
##### MD
# acoustic radiation_left
model2_15L <- lm(p25515_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_15L)
# acoustic radiation_right
model2_15R <- lm(p25516_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_15R)
# anterior thalamic radiation_left
model2_16L <- lm(p25517_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_16L)
# anterior thalamic radiation_right
model2_16R <- lm(p25518_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_16R)
# cingulate gyrus part of cingulum_left
model2_17L <- lm(p25519_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_17L)
# cingulate gyrus part of cingulum_right
model2_17R <- lm(p25520_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_17R)
# corticospinal tract_left
model2_18L <- lm(p25523_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_18L)
# corticospinal tract_right
model2_18R <- lm(p25524_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_18R)
# forceps major_left
model2_19L <- lm(p25525_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_19L)
# forceps major_right
model2_19R <- lm(p25526_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_19R)
# inferior fronto-occipital fasciculus_left
model2_20L <- lm(p25527_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_20L)
# inferior fronto-occipital fasciculus_right
model2_20R <- lm(p25528_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_20R)
# longitudinal fasciculus_left
model2_21L <- lm(p25529_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_21L)
# longitudinal fasciculus_right
model2_21R <- lm(p25530_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_21R)
# medial lemniscus_left
model2_22L <- lm(p25532_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_22L)
# medial lemniscus_right
model2_22R <- lm(p25533_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_22R)
# middle cerebellar peduncle
model2_23 <- lm(p25531_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_23)
# parahippocampal part of cingulum_left
model2_24L <- lm(p25521_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_24L)
# parahippocampal part of cingulum_right
model2_24R <- lm(p25522_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_24R)
# posterior thalamic radiation_left
model2_25L <- lm(p25534_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_25L)
# posterior thalamic radiation_right
model2_25R <- lm(p25535_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model2_25R)
# superior longitudinal fasciculus_left
model2_26L <- lm(p25536_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_26L)
# superior longitudinal fasciculus_right
model2_26R <- lm(p25537_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model2_26R)
# superior thalamic radiation_left
model2_27L <- lm(p25538_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_27L)
# superior thalamic radiation_right
model2_27R <- lm(p25539_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_27R)
# uncinate fasciculus_left
model2_28L <- lm(p25540_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_28L)
# uncinate fasciculus_right
model2_28R <- lm(p25541_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
summary(model2_28R)
# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model2_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  
  summary(model2_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_9)$coefficients["alcohol_pw", "Pr(>|t|)"], 
  summary(model2_10L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_10R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_11L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_11R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_12L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_12R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_13L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_13R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_14L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_14R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_15L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_15R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_16L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_16R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_17L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_17R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_18L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_18R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_19L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_19R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_20L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_20R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_21L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_21R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_22L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_22R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_23)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_24L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_24R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_25L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_25R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_26L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_26R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_27L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_27R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_28L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model2_28R)$coefficients["alcohol_pw", "Pr(>|t|)"]
)

# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)

##### 模型3
##### FA
# acoustic radiation_left
model3_1L <- lm(p25488_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_1L)
# acoustic radiation_right
model3_1R <- lm(p25489_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_1R)
# anterior thalamic radiation_left
model3_2L <- lm(p25490_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_2L)
# anterior thalamic radiation_right
model3_2R <- lm(p25491_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_2R)
# cingulate gyrus part of cingulum_left
model3_3L <- lm(p25492_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_3L)
# cingulate gyrus part of cingulum_right
model3_3R <- lm(p25493_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_3R)
# corticospinal tract_left
model3_4L <- lm(p25496_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_4L)
# corticospinal tract_right
model3_4R <- lm(p25497_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_4R)
# forceps major_left
model3_5L <- lm(p25498_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_5L)
# forceps major_right
model3_5R <- lm(p25499_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_5R)
# inferior fronto-occipital fasciculus_left
model3_6L <- lm(p25500_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_6L)
# inferior fronto-occipital fasciculus_right
model3_6R <- lm(p25501_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_6R)
# longitudinal fasciculus_left
model3_7L <- lm(p25502_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_7L)
# longitudinal fasciculus_right
model3_7R <- lm(p25503_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_7R)
# medial lemniscus_left
model3_8L <- lm(p25505_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_8L)
# medial lemniscus_right
model3_8R <- lm(p25506_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_8R)
# middle cerebellar peduncle
model3_9 <- lm(p25504_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_9)
# parahippocampal part of cingulum_left
model3_10L <- lm(p25494_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_10L)
# parahippocampal part of cingulum_right
model3_10R <- lm(p25495_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_10R)
# posterior thalamic radiation_left
model3_11L <- lm(p25507_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_11L)
# posterior thalamic radiation_right
model3_11R <- lm(p25508_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_11R)
# superior longitudinal fasciculus_left
model3_12L <- lm(p25509_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_12L)
# superior longitudinal fasciculus_right
model3_12R <- lm(p25510_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_12R)
# superior thalamic radiation_left
model3_13L <- lm(p25511_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_13L)
# superior thalamic radiation_right
model3_13R <- lm(p25512_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_13R)
# uncinate fasciculus_left
model3_14L <- lm(p25513_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_14L)
# uncinate fasciculus_right
model3_14R <- lm(p25514_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_14R)
##### 模型3
##### MD
# acoustic radiation_left
model3_15L <- lm(p25515_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_15L)
# acoustic radiation_right
model3_15R <- lm(p25516_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_15R)
# anterior thalamic radiation_left
model3_16L <- lm(p25517_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_16L)
# anterior thalamic radiation_right
model3_16R <- lm(p25518_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_16R)
# cingulate gyrus part of cingulum_left
model3_17L <- lm(p25519_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_17L)
# cingulate gyrus part of cingulum_right
model3_17R <- lm(p25520_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_17R)
# corticospinal tract_left
model3_18L <- lm(p25523_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_18L)
# corticospinal tract_right
model3_18R <- lm(p25524_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_18R)
# forceps major_left
model3_19L <- lm(p25525_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating +  triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_19L)
# forceps major_right
model3_19R <- lm(p25526_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_19R)
# inferior fronto-occipital fasciculus_left
model3_20L <- lm(p25527_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_20L)
# inferior fronto-occipital fasciculus_right
model3_20R <- lm(p25528_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_20R)
# longitudinal fasciculus_left
model3_21L <- lm(p25529_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_21L)
# longitudinal fasciculus_right
model3_21R <- lm(p25530_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_21R)
# medial lemniscus_left
model3_22L <- lm(p25532_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_22L)
# medial lemniscus_right
model3_22R <- lm(p25533_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_22R)
# middle cerebellar peduncle
model3_23 <- lm(p25531_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_23)
# parahippocampal part of cingulum_left
model3_24L <- lm(p25521_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_24L)
# parahippocampal part of cingulum_right
model3_24R <- lm(p25522_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_24R)
# posterior thalamic radiation_left
model3_25L <- lm(p25534_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_25L)
# posterior thalamic radiation_right
model3_25R <- lm(p25535_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_25R)
# superior longitudinal fasciculus_left
model3_26L <- lm(p25536_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality+ blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate , data = mri)
summary(model3_26L)
# superior longitudinal fasciculus_right
model3_26R <- lm(p25537_i2 ~ alcohol_pw + preference + with_food, data = mri)
summary(model3_26R)
# superior thalamic radiation_left
model3_27L <- lm(p25538_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_27L)
# superior thalamic radiation_right
model3_27R <- lm(p25539_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_27R)
# uncinate fasciculus_left
model3_28L <- lm(p25540_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_28L)
# uncinate fasciculus_right
model3_28R <- lm(p25541_i2 ~ alcohol_pw + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
summary(model3_28R)
# 提取每个模型的p值（以alcohol_pw为例）
p_values <- c(
  summary(model3_1L)$coefficients["alcohol_pw", "Pr(>|t|)"],  
  summary(model3_1R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_2L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_2R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_3L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_3R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_4L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_4R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_5L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_5R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_6L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_6R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_7L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_7R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_8L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_8R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_9)$coefficients["alcohol_pw", "Pr(>|t|)"], 
  summary(model3_10L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_10R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_11L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_11R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_12L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_12R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_13L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_13R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_14L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_14R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_15L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_15R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_16L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_16R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_17L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_17R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_18L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_18R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_19L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_19R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_20L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_20R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_21L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_21R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_22L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_22R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_23)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_24L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_24R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_25L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_25R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_26L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_26R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_27L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_27R)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_28L)$coefficients["alcohol_pw", "Pr(>|t|)"],
  summary(model3_28R)$coefficients["alcohol_pw", "Pr(>|t|)"]
)

# 进行FDR校正
fdr_adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(fdr_adjusted_p_values)