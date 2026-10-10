# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
getwd()
mri <- read.csv("mri.csv")
data <- read.csv("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\5-Cox\\data with prs3_r.csv")

# 加载必要的R包
library(survival)    # Cox回归分析
library(caret)
library(car)         # 计算VIF
library(carData)     # 计算VIF
library(broom)
library(dplyr)

##### 提取用于分析的MRI数据（与结局变量的Cox分析）
# 从 mri.csv 中删除那些不在 data.csv 中的 eid 对应的行
eids_in_data <- unique(data$eid)
mri <- mri %>% 
  filter(eid %in% eids_in_data)

##### 仅纳入第一次影像学数据进行分析
# 计数并删除MRI数据缺失的参与者人数
# 删除i3数据
i3_remove <- grep("i3", colnames(mri), ignore.case = TRUE)
mri <- mri[, -i3_remove]
is_all_blank <- function(row) {
  all(is.na(row) | row == "", na.rm = TRUE)
}
# 仅纳入皮质下体积、各脑区体积、FA、MD
# 删除第74列、第79-217列的数据
names(mri)
mri <- mri[, -c(74, 79:217)]
# 删除有任一mri数据缺失的人
rows_with_na <- rowSums(is.na(mri[, 2:131])) >= 1
num_rows_with_na <- sum(rows_with_na)
cat("有任一mri数据缺失的人数为：", num_rows_with_na, "\n")
mri <- mri[!rows_with_na, ]


# 将data的第5-26，30-34列的信息纳入mri信息中
# 从data.csv中提取需要合并的列
data_to_merge <- data %>%
  select(eid, 5:26, 30:34, 40)  # 提取eid、第5-26列
# 按照eid将data_to_merge合并到mri.csv中，只保留mri.csv中存在的eid
mri <- mri %>%
  left_join(data_to_merge, by = "eid")  # 使用左连接（left join）

write.csv(mri, "mri_mediation with time.csv", row.names = FALSE)


##### 读取用于分析的数据
# 设置并读取工作目录
setwd("D:\\读研倒计时\\课题\\UKB\\data\\merged data\\6-MRI")
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
mri$All_cause_stroke <- as.numeric(mri$All_cause_stroke)
mri$Hemorrhagic_stroke <- as.numeric(mri$Hemorrhagic_stroke)
mri$Ischemic_stroke <- as.numeric(mri$Ischemic_stroke)
mri$TIA <- as.numeric(mri$TIA)
str(mri)
mri$alcohol_pw <- scale(mri$alcohol_pw)

# 创建生存对象
mri$surv_obj <- Surv(mri$follow_up_time, mri$ACDs)

##### volume
##### 标准化处理
mri[, 2:77] <- scale(mri[, 2:77]) 

##### 模型1
##### 皮层下区域
# 检查模型1协变量之间的多重共线性
model_1 <- vif(lm(alcohol_pw ~ preference + with_food, data = mri))
print(model_1)
# 皮层下区域
model1_1L <- coxph(mri$surv_obj ~ alcohol_pw + p25023_i2 + preference + with_food, data = mri) # accumbens_left
model1_1R <- coxph(mri$surv_obj ~ alcohol_pw + p25024_i2 + preference + with_food, data = mri) # accumbens_right
model1_2L <- coxph(mri$surv_obj ~ alcohol_pw + p25021_i2 + preference + with_food, data = mri) # amygdala_left
model1_2R <- coxph(mri$surv_obj ~ alcohol_pw + p25022_i2 + preference + with_food, data = mri) # amygdala_right
model1_3L <- coxph(mri$surv_obj ~ alcohol_pw + p25013_i2 + preference + with_food, data = mri) # caudate_left
model1_3R <- coxph(mri$surv_obj ~ alcohol_pw + p25014_i2 + preference + with_food, data = mri) # caudate_right
model1_4L <- coxph(mri$surv_obj ~ alcohol_pw + p25019_i2 + preference + with_food, data = mri) # hippocampus_left
model1_4R <- coxph(mri$surv_obj ~ alcohol_pw + p25020_i2 + preference + with_food, data = mri) # hippocampus_right
model1_5L <- coxph(mri$surv_obj ~ alcohol_pw + p25017_i2 + preference + with_food, data = mri) # pallidum_left
model1_5R <- coxph(mri$surv_obj ~ alcohol_pw + p25018_i2 + preference + with_food, data = mri) # pallidum_right
model1_6L <- coxph(mri$surv_obj ~ alcohol_pw + p25015_i2 + preference + with_food, data = mri) # putamen_left
model1_6R <- coxph(mri$surv_obj ~ alcohol_pw + p25016_i2 + preference + with_food, data = mri) # putamen_right
model1_7L <- coxph(mri$surv_obj ~ alcohol_pw + p25011_i2 + preference + with_food, data = mri) # thalamus_left
model1_7R <- coxph(mri$surv_obj ~ alcohol_pw + p25012_i2 + preference + with_food, data = mri) # thalamus_right

# 各个脑区的体积
model1_8L <- coxph(mri$surv_obj ~ alcohol_pw + p27207_i2 + preference + with_food, data = mri) # cuneus_left
model1_8R <- coxph(mri$surv_obj ~ alcohol_pw + p27300_i2 + preference + with_food, data = mri) # cuneus_right
model1_9L <- coxph(mri$surv_obj ~ alcohol_pw + p27208_i2 + preference + with_food, data = mri) # entorhinal_left
model1_9R <- coxph(mri$surv_obj ~ alcohol_pw + p27301_i2 + preference + with_food, data = mri) # entorhinal_right
model1_10L <- coxph(mri$surv_obj ~ alcohol_pw + p27209_i2 + preference + with_food, data = mri) # fusiform_left
model1_10R <- coxph(mri$surv_obj ~ alcohol_pw + p27302_i2 + preference + with_food, data = mri) # fusiform_right
model1_11L <- coxph(mri$surv_obj ~ alcohol_pw + p27210_i2 + preference + with_food, data = mri) # inferiorparietal_left
model1_11R <- coxph(mri$surv_obj ~ alcohol_pw + p27303_i2 + preference + with_food, data = mri) # inferiorparietal_right
model1_12L <- coxph(mri$surv_obj ~ alcohol_pw + p27211_i2 + preference + with_food, data = mri) # inferiortemporal_left
model1_12R <- coxph(mri$surv_obj ~ alcohol_pw + p27304_i2 + preference + with_food, data = mri) # inferiortemporal_right
model1_13L <- coxph(mri$surv_obj ~ alcohol_pw + p27235_i2 + preference + with_food, data = mri) # insula_left
model1_13R <- coxph(mri$surv_obj ~ alcohol_pw + p27328_i2 + preference + with_food, data = mri) # insula_right
model1_14L <- coxph(mri$surv_obj ~ alcohol_pw + p27212_i2 + preference + with_food, data = mri) # isthmuscingulate_left
model1_14R <- coxph(mri$surv_obj ~ alcohol_pw + p27305_i2 + preference + with_food, data = mri) # isthmuscingulate_right
model1_15L <- coxph(mri$surv_obj ~ alcohol_pw + p27213_i2 + preference + with_food, data = mri) # lateraloccipital_left
model1_15R <- coxph(mri$surv_obj ~ alcohol_pw + p27306_i2 + preference + with_food, data = mri) # lateraloccipital_right
model1_16L <- coxph(mri$surv_obj ~ alcohol_pw + p27214_i2 + preference + with_food, data = mri) # lateralorbitofrontal_left
model1_16R <- coxph(mri$surv_obj ~ alcohol_pw + p27307_i2 + preference + with_food, data = mri) # lateralorbitofrontal_right
model1_17L <- coxph(mri$surv_obj ~ alcohol_pw + p27215_i2 + preference + with_food, data = mri) # lingual_left
model1_17R <- coxph(mri$surv_obj ~ alcohol_pw + p27308_i2 + preference + with_food, data = mri) # lingual_right
model1_18L <- coxph(mri$surv_obj ~ alcohol_pw + p27216_i2 + preference + with_food, data = mri) # medialorbitofrontal_left
model1_18R <- coxph(mri$surv_obj ~ alcohol_pw + p27309_i2 + preference + with_food, data = mri) # medialorbitofrontal_right
model1_19L <- coxph(mri$surv_obj ~ alcohol_pw + p27217_i2 + preference + with_food, data = mri) # middletemporal_left
model1_19R <- coxph(mri$surv_obj ~ alcohol_pw + p27310_i2 + preference + with_food, data = mri) # middletemporal_right
model1_20L <- coxph(mri$surv_obj ~ alcohol_pw + p27219_i2 + preference + with_food, data = mri) # paracentral_left
model1_20R <- coxph(mri$surv_obj ~ alcohol_pw + p27312_i2 + preference + with_food, data = mri) # paracentral_right
model1_21L <- coxph(mri$surv_obj ~ alcohol_pw + p27218_i2 + preference + with_food, data = mri) # parahippocampal_left
model1_21R <- coxph(mri$surv_obj ~ alcohol_pw + p27311_i2 + preference + with_food, data = mri) # parahippocampal_right
model1_22L <- coxph(mri$surv_obj ~ alcohol_pw + p27220_i2 + preference + with_food, data = mri) # parsopercularis_left
model1_22R <- coxph(mri$surv_obj ~ alcohol_pw + p27313_i2 + preference + with_food, data = mri) # parsopercularis_right
model1_23L <- coxph(mri$surv_obj ~ alcohol_pw + p27221_i2 + preference + with_food, data = mri) # parsorbitalis_left
model1_23R <- coxph(mri$surv_obj ~ alcohol_pw + p27314_i2 + preference + with_food, data = mri) # parsorbitalis_right
model1_24L <- coxph(mri$surv_obj ~ alcohol_pw + p27222_i2 + preference + with_food, data = mri) # parstriangularis_left
model1_24R <- coxph(mri$surv_obj ~ alcohol_pw + p27315_i2 + preference + with_food, data = mri) # parstriangularis_right
model1_25L <- coxph(mri$surv_obj ~ alcohol_pw + p27223_i2 + preference + with_food, data = mri) # pericalcarine_left
model1_25R <- coxph(mri$surv_obj ~ alcohol_pw + p27316_i2 + preference + with_food, data = mri) # pericalcarine_right
model1_26L <- coxph(mri$surv_obj ~ alcohol_pw + p27224_i2 + preference + with_food, data = mri) # postcentral_left
model1_26R <- coxph(mri$surv_obj ~ alcohol_pw + p27317_i2 + preference + with_food, data = mri) # postcentral_right
model1_27L <- coxph(mri$surv_obj ~ alcohol_pw + p27225_i2 + preference + with_food, data = mri) # posteriorcingulate_left
model1_27R <- coxph(mri$surv_obj ~ alcohol_pw + p27318_i2 + preference + with_food, data = mri) # posteriorcingulate_right
model1_28L <- coxph(mri$surv_obj ~ alcohol_pw + p27226_i2 + preference + with_food, data = mri) # precentral_left
model1_28R <- coxph(mri$surv_obj ~ alcohol_pw + p27319_i2 + preference + with_food, data = mri) # precentral_right
model1_29L <- coxph(mri$surv_obj ~ alcohol_pw + p27227_i2 + preference + with_food, data = mri) # precuneus_left
model1_29R <- coxph(mri$surv_obj ~ alcohol_pw + p27320_i2 + preference + with_food, data = mri) # precuneus_right
model1_30L <- coxph(mri$surv_obj ~ alcohol_pw + p27228_i2 + preference + with_food, data = mri) # rostralanteriorcingulate_left
model1_30R <- coxph(mri$surv_obj ~ alcohol_pw + p27321_i2 + preference + with_food, data = mri) # rostralanteriorcingulate_right
model1_31L <- coxph(mri$surv_obj ~ alcohol_pw + p27229_i2 + preference + with_food, data = mri) # rostralmiddlefrontal_left
model1_31R <- coxph(mri$surv_obj ~ alcohol_pw + p27322_i2 + preference + with_food, data = mri) # rostralmiddlefrontal_right
model1_32L <- coxph(mri$surv_obj ~ alcohol_pw + p27230_i2 + preference + with_food, data = mri) # superiorfrontal_left
model1_32R <- coxph(mri$surv_obj ~ alcohol_pw + p27323_i2 + preference + with_food, data = mri) # superiorfrontal_right
model1_33L <- coxph(mri$surv_obj ~ alcohol_pw + p27231_i2 + preference + with_food, data = mri) # superiorparietal_left
model1_33R <- coxph(mri$surv_obj ~ alcohol_pw + p27324_i2 + preference + with_food, data = mri) # superiorparietal_right
model1_34L <- coxph(mri$surv_obj ~ alcohol_pw + p27232_i2 + preference + with_food, data = mri) # superiortemporal_left
model1_34R <- coxph(mri$surv_obj ~ alcohol_pw + p27325_i2 + preference + with_food, data = mri) # superiortemporal_right
model1_35L <- coxph(mri$surv_obj ~ alcohol_pw + p27233_i2 + preference + with_food, data = mri) # supramarginal_left
model1_35R <- coxph(mri$surv_obj ~ alcohol_pw + p27326_i2 + preference + with_food, data = mri) # supramarginal_right
model1_36L <- coxph(mri$surv_obj ~ alcohol_pw + p27234_i2 + preference + with_food, data = mri) # transversetemporal_left
model1_36R <- coxph(mri$surv_obj ~ alcohol_pw + p27327_i2 + preference + with_food, data = mri) # transversetemporal_right
model1_37L <- coxph(mri$surv_obj ~ alcohol_pw + p27205_i2 + preference + with_food, data = mri) # caudalanteriorcingulate_left
model1_37R <- coxph(mri$surv_obj ~ alcohol_pw + p27298_i2 + preference + with_food, data = mri) # caudalanteriorcingulate_right
model1_38L <- coxph(mri$surv_obj ~ alcohol_pw + p27206_i2 + preference + with_food, data = mri) # caudalmiddlefrontal_left
model1_38R <- coxph(mri$surv_obj ~ alcohol_pw + p27299_i2 + preference + with_food, data = mri) # caudalmiddlefrontal_right

# 输出结果
summary(model1_1L)
summary(model1_1R)
summary(model1_2L)
summary(model1_2R)
summary(model1_3L)
summary(model1_3R)
summary(model1_4L)
summary(model1_4R)
summary(model1_5L)
summary(model1_5R)
summary(model1_6L)
summary(model1_6R)
summary(model1_7L)
summary(model1_7R)
summary(model1_8L)
summary(model1_8R)
summary(model1_9L)
summary(model1_9R)
summary(model1_10L)
summary(model1_10R)
summary(model1_11L)
summary(model1_11R)
summary(model1_12L)
summary(model1_12R)
summary(model1_13L)
summary(model1_13R)
summary(model1_14L)
summary(model1_14R)
summary(model1_15L)
summary(model1_15R)
summary(model1_16L)
summary(model1_16R)
summary(model1_17L)
summary(model1_17R)
summary(model1_18L)
summary(model1_18R)
summary(model1_19L)
summary(model1_19R)
summary(model1_20L)
summary(model1_20R)
summary(model1_21L)
summary(model1_21R)
summary(model1_22L)
summary(model1_22R)
summary(model1_23L)
summary(model1_23R)
summary(model1_24L)
summary(model1_24R)
summary(model1_25L)
summary(model1_25R)
summary(model1_26L)
summary(model1_26R)
summary(model1_27L)
summary(model1_27R)
summary(model1_28L)
summary(model1_28R)
summary(model1_29L)
summary(model1_29R)
summary(model1_30L)
summary(model1_30R)
summary(model1_31L)
summary(model1_31R)
summary(model1_32L)
summary(model1_32R)
summary(model1_33L)
summary(model1_33R)
summary(model1_34L)
summary(model1_34R)
summary(model1_35L)
summary(model1_35R)
summary(model1_36L)
summary(model1_36R)
summary(model1_37L)
summary(model1_37R)
summary(model1_38L)
summary(model1_38R)


##### 模型2
# 检查模型2协变量的多重线性关系
model_2 <- vif(lm(alcohol_pw ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri))
print(model_2)
# 基底节区域
model2_1L <- coxph(mri$surv_obj ~ alcohol_pw + p25023_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # accumbens_left
model2_1R <- coxph(mri$surv_obj ~ alcohol_pw + p25024_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # accumbens_right
model2_2L <- coxph(mri$surv_obj ~ alcohol_pw + p25021_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # amygdala_left
model2_2R <- coxph(mri$surv_obj ~ alcohol_pw + p25022_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # amygdala_right
model2_3L <- coxph(mri$surv_obj ~ alcohol_pw + p25013_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudate_left
model2_3R <- coxph(mri$surv_obj ~ alcohol_pw + p25014_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudate_right
model2_4L <- coxph(mri$surv_obj ~ alcohol_pw + p25019_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # hippocampus_left
model2_4R <- coxph(mri$surv_obj ~ alcohol_pw + p25020_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # hippocampus_right
model2_5L <- coxph(mri$surv_obj ~ alcohol_pw + p25017_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # pallidum_left
model2_5R <- coxph(mri$surv_obj ~ alcohol_pw + p25018_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # pallidum_right
model2_6L <- coxph(mri$surv_obj ~ alcohol_pw + p25015_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # putamen_left
model2_6R <- coxph(mri$surv_obj ~ alcohol_pw + p25016_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # putamen_right
model2_7L <- coxph(mri$surv_obj ~ alcohol_pw + p25011_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # thalamus_left
model2_7R <- coxph(mri$surv_obj ~ alcohol_pw + p25012_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # thalamus_right

# 各个脑区的体积
model2_8L <- coxph(mri$surv_obj ~ alcohol_pw + p27207_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # cuneus_left
model2_8R <- coxph(mri$surv_obj ~ alcohol_pw + p27300_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # cuneus_right
model2_9L <- coxph(mri$surv_obj ~ alcohol_pw + p27208_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # entorhinal_left
model2_9R <- coxph(mri$surv_obj ~ alcohol_pw + p27301_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # entorhinal_right
model2_10L <- coxph(mri$surv_obj ~ alcohol_pw + p27209_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # fusiform_left
model2_10R <- coxph(mri$surv_obj ~ alcohol_pw + p27302_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # fusiform_right
model2_11L <- coxph(mri$surv_obj ~ alcohol_pw + p27210_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # inferiorparietal_left
model2_11R <- coxph(mri$surv_obj ~ alcohol_pw + p27303_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # inferiorparietal_right
model2_12L <- coxph(mri$surv_obj ~ alcohol_pw + p27211_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # inferiortemporal_left
model2_12R <- coxph(mri$surv_obj ~ alcohol_pw + p27304_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # inferiortemporal_right
model2_13L <- coxph(mri$surv_obj ~ alcohol_pw + p27235_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # insula_left
model2_13R <- coxph(mri$surv_obj ~ alcohol_pw + p27328_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # insula_right
model2_14L <- coxph(mri$surv_obj ~ alcohol_pw + p27212_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # isthmuscingulate_left
model2_14R <- coxph(mri$surv_obj ~ alcohol_pw + p27305_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # isthmuscingulate_right
model2_15L <- coxph(mri$surv_obj ~ alcohol_pw + p27213_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lateraloccipital_left
model2_15R <- coxph(mri$surv_obj ~ alcohol_pw + p27306_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lateraloccipital_right
model2_16L <- coxph(mri$surv_obj ~ alcohol_pw + p27214_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lateralorbitofrontal_left
model2_16R <- coxph(mri$surv_obj ~ alcohol_pw + p27307_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lateralorbitofrontal_right
model2_17L <- coxph(mri$surv_obj ~ alcohol_pw + p27215_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lingual_left
model2_17R <- coxph(mri$surv_obj ~ alcohol_pw + p27308_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # lingual_right
model2_18L <- coxph(mri$surv_obj ~ alcohol_pw + p27216_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # medialorbitofrontal_left
model2_18R <- coxph(mri$surv_obj ~ alcohol_pw + p27309_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # medialorbitofrontal_right
model2_19L <- coxph(mri$surv_obj ~ alcohol_pw + p27217_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # middletemporal_left
model2_19R <- coxph(mri$surv_obj ~ alcohol_pw + p27310_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # middletemporal_right
model2_20L <- coxph(mri$surv_obj ~ alcohol_pw + p27219_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # paracentral_left
model2_20R <- coxph(mri$surv_obj ~ alcohol_pw + p27312_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # paracentral_right
model2_21L <- coxph(mri$surv_obj ~ alcohol_pw + p27218_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parahippocampal_left
model2_21R <- coxph(mri$surv_obj ~ alcohol_pw + p27311_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parahippocampal_right
model2_22L <- coxph(mri$surv_obj ~ alcohol_pw + p27220_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parsopercularis_left
model2_22R <- coxph(mri$surv_obj ~ alcohol_pw + p27313_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parsopercularis_right
model2_23L <- coxph(mri$surv_obj ~ alcohol_pw + p27221_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parsorbitalis_left
model2_23R <- coxph(mri$surv_obj ~ alcohol_pw + p27314_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parsorbitalis_right
model2_24L <- coxph(mri$surv_obj ~ alcohol_pw + p27222_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parstriangularis_left
model2_24R <- coxph(mri$surv_obj ~ alcohol_pw + p27315_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # parstriangularis_right
model2_25L <- coxph(mri$surv_obj ~ alcohol_pw + p27223_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # pericalcarine_left
model2_25R <- coxph(mri$surv_obj ~ alcohol_pw + p27316_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # pericalcarine_right
model2_26L <- coxph(mri$surv_obj ~ alcohol_pw + p27224_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # postcentral_left
model2_26R <- coxph(mri$surv_obj ~ alcohol_pw + p27317_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # postcentral_right
model2_27L <- coxph(mri$surv_obj ~ alcohol_pw + p27225_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # posteriorcingulate_left
model2_27R <- coxph(mri$surv_obj ~ alcohol_pw + p27318_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # posteriorcingulate_right
model2_28L <- coxph(mri$surv_obj ~ alcohol_pw + p27226_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # precentral_left
model2_28R <- coxph(mri$surv_obj ~ alcohol_pw + p27319_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # precentral_right
model2_29L <- coxph(mri$surv_obj ~ alcohol_pw + p27227_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # precuneus_left
model2_29R <- coxph(mri$surv_obj ~ alcohol_pw + p27320_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # precuneus_right
model2_30L <- coxph(mri$surv_obj ~ alcohol_pw + p27228_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # rostralanteriorcingulate_left
model2_30R <- coxph(mri$surv_obj ~ alcohol_pw + p27321_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # rostralanteriorcingulate_right
model2_31L <- coxph(mri$surv_obj ~ alcohol_pw + p27229_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # rostralmiddlefrontal_left
model2_31R <- coxph(mri$surv_obj ~ alcohol_pw + p27322_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # rostralmiddlefrontal_right
model2_32L <- coxph(mri$surv_obj ~ alcohol_pw + p27230_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiorfrontal_left
model2_32R <- coxph(mri$surv_obj ~ alcohol_pw + p27323_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiorfrontal_right
model2_33L <- coxph(mri$surv_obj ~ alcohol_pw + p27231_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiorparietal_left
model2_33R <- coxph(mri$surv_obj ~ alcohol_pw + p27324_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiorparietal_right
model2_34L <- coxph(mri$surv_obj ~ alcohol_pw + p27232_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiortemporal_left
model2_34R <- coxph(mri$surv_obj ~ alcohol_pw + p27325_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # superiortemporal_right
model2_35L <- coxph(mri$surv_obj ~ alcohol_pw + p27233_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # supramarginal_left
model2_35R <- coxph(mri$surv_obj ~ alcohol_pw + p27326_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # supramarginal_right
model2_36L <- coxph(mri$surv_obj ~ alcohol_pw + p27234_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # transversetemporal_left
model2_36R <- coxph(mri$surv_obj ~ alcohol_pw + p27327_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # transversetemporal_right
model2_37L <- coxph(mri$surv_obj ~ alcohol_pw + p27205_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudalanteriorcingulate_left
model2_37R <- coxph(mri$surv_obj ~ alcohol_pw + p27298_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudalanteriorcingulate_right
model2_38L <- coxph(mri$surv_obj ~ alcohol_pw + p27206_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudalmiddlefrontal_left
model2_38R <- coxph(mri$surv_obj ~ alcohol_pw + p27299_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # caudalmiddlefrontal_right


# 输出结果
summary(model2_1L)
summary(model2_1R)
summary(model2_2L)
summary(model2_2R)
summary(model2_3L)
summary(model2_3R)
summary(model2_4L)
summary(model2_4R)
summary(model2_5L)
summary(model2_5R)
summary(model2_6L)
summary(model2_6R)
summary(model2_7L)
summary(model2_7R)
summary(model2_8L)
summary(model2_8R)
summary(model2_9L)
summary(model2_9R)
summary(model2_10L)
summary(model2_10R)
summary(model2_11L)
summary(model2_11R)
summary(model2_12L)
summary(model2_12R)
summary(model2_13L)
summary(model2_13R)
summary(model2_14L)
summary(model2_14R)
summary(model2_15L)
summary(model2_15R)
summary(model2_16L)
summary(model2_16R)
summary(model2_17L)
summary(model2_17R)
summary(model2_18L)
summary(model2_18R)
summary(model2_19L)
summary(model2_19R)
summary(model2_20L)
summary(model2_20R)
summary(model2_21L)
summary(model2_21R)
summary(model2_22L)
summary(model2_22R)
summary(model2_23L)
summary(model2_23R)
summary(model2_24L)
summary(model2_24R)
summary(model2_25L)
summary(model2_25R)
summary(model2_26L)
summary(model2_26R)
summary(model2_27L)
summary(model2_27R)
summary(model2_28L)
summary(model2_28R)
summary(model2_29L)
summary(model2_29R)
summary(model2_30L)
summary(model2_30R)
summary(model2_31L)
summary(model2_31R)
summary(model2_32L)
summary(model2_32R)
summary(model2_33L)
summary(model2_33R)
summary(model2_34L)
summary(model2_34R)
summary(model2_35L)
summary(model2_35R)
summary(model2_36L)
summary(model2_36R)
summary(model2_37L)
summary(model2_37R)
summary(model2_38L)
summary(model2_38R)


##### 模型3
##### 皮层下区域
# 检查模型3协变量之间的多重共线性
model_3 <- vif(lm(alcohol_pw ~ preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri))
print(model_3)
# 皮层下区域
model3_1L <- coxph(mri$surv_obj ~ alcohol_pw + p25023_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # accumbens_left
model3_1R <- coxph(mri$surv_obj ~ alcohol_pw + p25024_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # accumbens_right
model3_2L <- coxph(mri$surv_obj ~ alcohol_pw + p25021_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # amygdala_left
model3_2R <- coxph(mri$surv_obj ~ alcohol_pw + p25022_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # amygdala_right
model3_3L <- coxph(mri$surv_obj ~ alcohol_pw + p25013_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudate_left
model3_3R <- coxph(mri$surv_obj ~ alcohol_pw + p25014_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudate_right
model3_4L <- coxph(mri$surv_obj ~ alcohol_pw + p25019_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # hippocampus_left
model3_4R <- coxph(mri$surv_obj ~ alcohol_pw + p25020_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # hippocampus_right
model3_5L <- coxph(mri$surv_obj ~ alcohol_pw + p25017_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # pallidum_left
model3_5R <- coxph(mri$surv_obj ~ alcohol_pw + p25018_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # pallidum_right
model3_6L <- coxph(mri$surv_obj ~ alcohol_pw + p25015_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # putamen_left
model3_6R <- coxph(mri$surv_obj ~ alcohol_pw + p25016_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # putamen_right
model3_7L <- coxph(mri$surv_obj ~ alcohol_pw + p25011_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # thalamus_left
model3_7R <- coxph(mri$surv_obj ~ alcohol_pw + p25012_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # thalamus_right

# 各脑区的体积
model3_8L <- coxph(mri$surv_obj ~ alcohol_pw + p27207_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # cuneus_left
model3_8R <- coxph(mri$surv_obj ~ alcohol_pw + p27300_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # cuneus_right
model3_9L <- coxph(mri$surv_obj ~ alcohol_pw + p27208_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # entorhinal_left
model3_9R <- coxph(mri$surv_obj ~ alcohol_pw + p27301_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # entorhinal_right
model3_10L <- coxph(mri$surv_obj ~ alcohol_pw + p27209_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # fusiform_left
model3_10R <- coxph(mri$surv_obj ~ alcohol_pw + p27302_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # fusiform_right
model3_11L <- coxph(mri$surv_obj ~ alcohol_pw + p27210_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # inferiorparietal_left
model3_11R <- coxph(mri$surv_obj ~ alcohol_pw + p27303_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # inferiorparietal_right
model3_12L <- coxph(mri$surv_obj ~ alcohol_pw + p27211_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # inferiortemporal_left
model3_12R <- coxph(mri$surv_obj ~ alcohol_pw + p27304_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # inferiortemporal_right
model3_13L <- coxph(mri$surv_obj ~ alcohol_pw + p27235_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # insula_left
model3_13R <- coxph(mri$surv_obj ~ alcohol_pw + p27328_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # insula_right
model3_14L <- coxph(mri$surv_obj ~ alcohol_pw + p27212_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # isthmuscingulate_left
model3_14R <- coxph(mri$surv_obj ~ alcohol_pw + p27305_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # isthmuscingulate_right
model3_15L <- coxph(mri$surv_obj ~ alcohol_pw + p27213_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lateraloccipital_left
model3_15R <- coxph(mri$surv_obj ~ alcohol_pw + p27306_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lateraloccipital_right
model3_16L <- coxph(mri$surv_obj ~ alcohol_pw + p27214_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lateralorbitofrontal_left
model3_16R <- coxph(mri$surv_obj ~ alcohol_pw + p27307_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lateralorbitofrontal_right
model3_17L <- coxph(mri$surv_obj ~ alcohol_pw + p27215_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lingual_left
model3_17R <- coxph(mri$surv_obj ~ alcohol_pw + p27308_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # lingual_right
model3_18L <- coxph(mri$surv_obj ~ alcohol_pw + p27216_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # medialorbitofrontal_left
model3_18R <- coxph(mri$surv_obj ~ alcohol_pw + p27309_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # medialorbitofrontal_right
model3_19L <- coxph(mri$surv_obj ~ alcohol_pw + p27217_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # middletemporal_left
model3_19R <- coxph(mri$surv_obj ~ alcohol_pw + p27310_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # middletemporal_right
model3_20L <- coxph(mri$surv_obj ~ alcohol_pw + p27219_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # paracentral_left
model3_20R <- coxph(mri$surv_obj ~ alcohol_pw + p27312_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # paracentral_right
model3_21L <- coxph(mri$surv_obj ~ alcohol_pw + p27218_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parahippocampal_left
model3_21R <- coxph(mri$surv_obj ~ alcohol_pw + p27311_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parahippocampal_right
model3_22L <- coxph(mri$surv_obj ~ alcohol_pw + p27220_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parsopercularis_left
model3_22R <- coxph(mri$surv_obj ~ alcohol_pw + p27313_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parsopercularis_right
model3_23L <- coxph(mri$surv_obj ~ alcohol_pw + p27221_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parsorbitalis_left
model3_23R <- coxph(mri$surv_obj ~ alcohol_pw + p27314_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parsorbitalis_right
model3_24L <- coxph(mri$surv_obj ~ alcohol_pw + p27222_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parstriangularis_left
model3_24R <- coxph(mri$surv_obj ~ alcohol_pw + p27315_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # parstriangularis_right
model3_25L <- coxph(mri$surv_obj ~ alcohol_pw + p27223_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # pericalcarine_left
model3_25R <- coxph(mri$surv_obj ~ alcohol_pw + p27316_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # pericalcarine_right
model3_26L <- coxph(mri$surv_obj ~ alcohol_pw + p27224_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # postcentral_left
model3_26R <- coxph(mri$surv_obj ~ alcohol_pw + p27317_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # postcentral_right
model3_27L <- coxph(mri$surv_obj ~ alcohol_pw + p27225_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # posteriorcingulate_left
model3_27R <- coxph(mri$surv_obj ~ alcohol_pw + p27318_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # posteriorcingulate_right
model3_28L <- coxph(mri$surv_obj ~ alcohol_pw + p27226_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # precentral_left
model3_28R <- coxph(mri$surv_obj ~ alcohol_pw + p27319_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # precentral_right
model3_29L <- coxph(mri$surv_obj ~ alcohol_pw + p27227_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # precuneus_left
model3_29R <- coxph(mri$surv_obj ~ alcohol_pw + p27320_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # precuneus_right
model3_30L <- coxph(mri$surv_obj ~ alcohol_pw + p27228_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # rostralanteriorcingulate_left
model3_30R <- coxph(mri$surv_obj ~ alcohol_pw + p27321_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # rostralanteriorcingulate_right
model3_31L <- coxph(mri$surv_obj ~ alcohol_pw + p27229_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # rostralmiddlefrontal_left
model3_31R <- coxph(mri$surv_obj ~ alcohol_pw + p27322_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # rostralmiddlefrontal_right
model3_32L <- coxph(mri$surv_obj ~ alcohol_pw + p27230_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiorfrontal_left
model3_32R <- coxph(mri$surv_obj ~ alcohol_pw + p27323_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiorfrontal_right
model3_33L <- coxph(mri$surv_obj ~ alcohol_pw + p27231_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiorparietal_left
model3_33R <- coxph(mri$surv_obj ~ alcohol_pw + p27324_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiorparietal_right
model3_34L <- coxph(mri$surv_obj ~ alcohol_pw + p27232_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiortemporal_left
model3_34R <- coxph(mri$surv_obj ~ alcohol_pw + p27325_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # superiortemporal_right
model3_35L <- coxph(mri$surv_obj ~ alcohol_pw + p27233_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # supramarginal_left
model3_35R <- coxph(mri$surv_obj ~ alcohol_pw + p27326_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # supramarginal_right
model3_36L <- coxph(mri$surv_obj ~ alcohol_pw + p27234_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # transversetemporal_left
model3_36R <- coxph(mri$surv_obj ~ alcohol_pw + p27327_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # transversetemporal_right
model3_37L <- coxph(mri$surv_obj ~ alcohol_pw + p27205_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudalanteriorcingulate_left
model3_37R <- coxph(mri$surv_obj ~ alcohol_pw + p27298_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudalanteriorcingulate_right
model3_38L <- coxph(mri$surv_obj ~ alcohol_pw + p27206_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudalmiddlefrontal_left
model3_38R <- coxph(mri$surv_obj ~ alcohol_pw + p27299_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # caudalmiddlefrontal_right


# 输出结果
summary(model3_1L)
summary(model3_1R)
summary(model3_2L)
summary(model3_2R)
summary(model3_3L)
summary(model3_3R)
summary(model3_4L)
summary(model3_4R)
summary(model3_5L)
summary(model3_5R)
summary(model3_6L)
summary(model3_6R)
summary(model3_7L)
summary(model3_7R)
summary(model3_8L)
summary(model3_8R)
summary(model3_9L)
summary(model3_9R)
summary(model3_10L)
summary(model3_10R)
summary(model3_11L)
summary(model3_11R)
summary(model3_12L)
summary(model3_12R)
summary(model3_13L)
summary(model3_13R)
summary(model3_14L)
summary(model3_14R)
summary(model3_15L)
summary(model3_15R)
summary(model3_16L)
summary(model3_16R)
summary(model3_17L)
summary(model3_17R)
summary(model3_18L)
summary(model3_18R)
summary(model3_19L)
summary(model3_19R)
summary(model3_20L)
summary(model3_20R)
summary(model3_21L)
summary(model3_21R)
summary(model3_22L)
summary(model3_22R)
summary(model3_23L)
summary(model3_23R)
summary(model3_24L)
summary(model3_24R)
summary(model3_25L)
summary(model3_25R)
summary(model3_26L)
summary(model3_26R)
summary(model3_27L)
summary(model3_27R)
summary(model3_28L)
summary(model3_28R)
summary(model3_29L)
summary(model3_29R)
summary(model3_30L)
summary(model3_30R)
summary(model3_31L)
summary(model3_31R)
summary(model3_32L)
summary(model3_32R)
summary(model3_33L)
summary(model3_33R)
summary(model3_34L)
summary(model3_34R)
summary(model3_35L)
summary(model3_35R)
summary(model3_36L)
summary(model3_36R)
summary(model3_37L)
summary(model3_37R)
summary(model3_38L)
summary(model3_38R)



##### FA MD 
# 创建生存对象
mri$surv_obj <- Surv(time = mri$time, event = mri$status)
#####FA、MD标准化处理
mri[, 78:104] <- scale(mri[, 78:104])  # 第一组特征标准化
mri[, 105:131] <- scale(mri[, 105:131])  # 第二组特征标准化


### 模型1 
# 白质纤维束FA指标分析
# acoustic radiation
model1_1L <- coxph(surv_obj ~ alcohol_pw + p25488_i2 + preference + with_food, data = mri) # left
model1_1R <- coxph(surv_obj ~ alcohol_pw + p25489_i2 + preference + with_food, data = mri) # right
# anterior thalamic radiation
model1_2L <- coxph(surv_obj ~ alcohol_pw + p25490_i2 + preference + with_food, data = mri) # left
model1_2R <- coxph(surv_obj ~ alcohol_pw + p25491_i2 + preference + with_food, data = mri) # right
# cingulate gyrus part of cingulum
model1_3L <- coxph(surv_obj ~ alcohol_pw + p25492_i2 + preference + with_food, data = mri) # left
model1_3R <- coxph(surv_obj ~ alcohol_pw + p25493_i2 + preference + with_food, data = mri) # right
# corticospinal tract
model1_4L <- coxph(surv_obj ~ alcohol_pw + p25496_i2 + preference + with_food, data = mri) # left
model1_4R <- coxph(surv_obj ~ alcohol_pw + p25497_i2 + preference + with_food, data = mri) # right
# forceps major
model1_5L <- coxph(surv_obj ~ alcohol_pw + p25498_i2 + preference + with_food, data = mri) # left
model1_5R <- coxph(surv_obj ~ alcohol_pw + p25499_i2 + preference + with_food, data = mri) # right
# inferior fronto-occipital fasciculus
model1_6L <- coxph(surv_obj ~ alcohol_pw + p25500_i2 + preference + with_food, data = mri) # left
model1_6R <- coxph(surv_obj ~ alcohol_pw + p25501_i2 + preference + with_food, data = mri) # right
# longitudinal fasciculus
model1_7L <- coxph(surv_obj ~ alcohol_pw + p25502_i2 + preference + with_food, data = mri) # left
model1_7R <- coxph(surv_obj ~ alcohol_pw + p25503_i2 + preference + with_food, data = mri) # right
# medial lemniscus
model1_8L <- coxph(surv_obj ~ alcohol_pw + p25505_i2 + preference + with_food, data = mri) # left
model1_8R <- coxph(surv_obj ~ alcohol_pw + p25506_i2 + preference + with_food, data = mri) # right
# middle cerebellar peduncle
model1_9 <- coxph(surv_obj ~ alcohol_pw + p25504_i2 + preference + with_food, data = mri)
# parahippocampal part of cingulum
model1_10L <- coxph(surv_obj ~ alcohol_pw + p25494_i2 + preference + with_food, data = mri) # left
model1_10R <- coxph(surv_obj ~ alcohol_pw + p25495_i2 + preference + with_food, data = mri) # right
# posterior thalamic radiation
model1_11L <- coxph(surv_obj ~ alcohol_pw + p25507_i2 + preference + with_food, data = mri) # left
model1_11R <- coxph(surv_obj ~ alcohol_pw + p25508_i2 + preference + with_food, data = mri) # right
# superior longitudinal fasciculus
model1_12L <- coxph(surv_obj ~ alcohol_pw + p25509_i2 + preference + with_food, data = mri) # left
model1_12R <- coxph(surv_obj ~ alcohol_pw + p25510_i2 + preference + with_food, data = mri) # right
# superior thalamic radiation
model1_13L <- coxph(surv_obj ~ alcohol_pw + p25511_i2 + preference + with_food, data = mri) # left
model1_13R <- coxph(surv_obj ~ alcohol_pw + p25512_i2 + preference + with_food, data = mri) # right
# uncinate fasciculus
model1_14L <- coxph(surv_obj ~ alcohol_pw + p25513_i2 + preference + with_food, data = mri) # left
model1_14R <- coxph(surv_obj ~ alcohol_pw + p25514_i2 + preference + with_food, data = mri) # right
# 白质纤维束MD指标分析
# acoustic radiation
model1_15L <- coxph(surv_obj ~ alcohol_pw + p25515_i2 + preference + with_food, data = mri) # left
model1_15R <- coxph(surv_obj ~ alcohol_pw + p25516_i2 + preference + with_food, data = mri) # right
# anterior thalamic radiation
model1_16L <- coxph(surv_obj ~ alcohol_pw + p25517_i2 + preference + with_food, data = mri) # left
model1_16R <- coxph(surv_obj ~ alcohol_pw + p25518_i2 + preference + with_food, data = mri) # right
# cingulate gyrus part of cingulum
model1_17L <- coxph(surv_obj ~ alcohol_pw + p25519_i2 + preference + with_food, data = mri) # left
model1_17R <- coxph(surv_obj ~ alcohol_pw + p25520_i2 + preference + with_food, data = mri) # right
# corticospinal tract
model1_18L <- coxph(surv_obj ~ alcohol_pw + p25523_i2 + preference + with_food, data = mri) # left
model1_18R <- coxph(surv_obj ~ alcohol_pw + p25524_i2 + preference + with_food, data = mri) # right
# forceps major
model1_19L <- coxph(surv_obj ~ alcohol_pw + p25525_i2 + preference + with_food, data = mri) # left
model1_19R <- coxph(surv_obj ~ alcohol_pw + p25526_i2 + preference + with_food, data = mri) # right
# inferior fronto-occipital fasciculus
model1_20L <- coxph(surv_obj ~ alcohol_pw + p25527_i2 + preference + with_food, data = mri) # left
model1_20R <- coxph(surv_obj ~ alcohol_pw + p25528_i2 + preference + with_food, data = mri) # right
# longitudinal fasciculus
model1_21L <- coxph(surv_obj ~ alcohol_pw + p25529_i2 + preference + with_food, data = mri) # left
model1_21R <- coxph(surv_obj ~ alcohol_pw + p25530_i2 + preference + with_food, data = mri) # right
# medial lemniscus
model1_22L <- coxph(surv_obj ~ alcohol_pw + p25532_i2 + preference + with_food, data = mri) # left
model1_22R <- coxph(surv_obj ~ alcohol_pw + p25533_i2 + preference + with_food, data = mri) # right
# middle cerebellar peduncle
model1_23 <- coxph(surv_obj ~ alcohol_pw + p25531_i2 + preference + with_food, data = mri)
# parahippocampal part of cingulum
model1_24L <- coxph(surv_obj ~ alcohol_pw + p25521_i2 + preference + with_food, data = mri) # left
model1_24R <- coxph(surv_obj ~ alcohol_pw + p25522_i2 + preference + with_food, data = mri) # right
# posterior thalamic radiation
model1_25L <- coxph(surv_obj ~ alcohol_pw + p25534_i2 + preference + with_food, data = mri) # left
model1_25R <- coxph(surv_obj ~ alcohol_pw + p25535_i2 + preference + with_food, data = mri) # right
# superior longitudinal fasciculus
model1_26L <- coxph(surv_obj ~ alcohol_pw + p25536_i2 + preference + with_food, data = mri) # left
model1_26R <- coxph(surv_obj ~ alcohol_pw + p25537_i2 + preference + with_food, data = mri) # right
# superior thalamic radiation
model1_27L <- coxph(surv_obj ~ alcohol_pw + p25538_i2 + preference + with_food, data = mri) # left
model1_27R <- coxph(surv_obj ~ alcohol_pw + p25539_i2 + preference + with_food, data = mri) # right
# uncinate fasciculus
model1_28L <- coxph(surv_obj ~ alcohol_pw + p25540_i2 + preference + with_food, data = mri) # left
model1_28R <- coxph(surv_obj ~ alcohol_pw + p25541_i2 + preference + with_food, data = mri) # right

# 模型数据归纳
summary(model1_1L)
summary(model1_1R)
summary(model1_2L)
summary(model1_2R)
summary(model1_3L)
summary(model1_3R)
summary(model1_4L)
summary(model1_4R)
summary(model1_5L)
summary(model1_5R)
summary(model1_6L)
summary(model1_6R)
summary(model1_7L)
summary(model1_7R)
summary(model1_8L)
summary(model1_8R)
summary(model1_9)
summary(model1_10L)
summary(model1_10R)
summary(model1_11L)
summary(model1_11R)
summary(model1_12L)
summary(model1_12R)
summary(model1_13L)
summary(model1_13R)
summary(model1_14L)
summary(model1_14R)
summary(model1_15L)
summary(model1_15R)
summary(model1_16L)
summary(model1_16R)
summary(model1_17L)
summary(model1_17R)
summary(model1_18L)
summary(model1_18R)
summary(model1_19L)
summary(model1_19R)
summary(model1_20L)
summary(model1_20R)
summary(model1_21L)
summary(model1_21R)
summary(model1_22L)
summary(model1_22R)
summary(model1_23)
summary(model1_24L)
summary(model1_24R)
summary(model1_25L)
summary(model1_25R)
summary(model1_26L)
summary(model1_26R)
summary(model1_27L)
summary(model1_27R)
summary(model1_28L)
summary(model1_28R)


### 模型2 
# 白质纤维束FA指标分析
# acoustic radiation
model2_1L <- coxph(surv_obj ~ alcohol_pw + p25488_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_1R <- coxph(surv_obj ~ alcohol_pw + p25489_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# anterior thalamic radiation
model2_2L <- coxph(surv_obj ~ alcohol_pw + p25490_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_2R <- coxph(surv_obj ~ alcohol_pw + p25491_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# cingulate gyrus part of cingulum
model2_3L <- coxph(surv_obj ~ alcohol_pw + p25492_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_3R <- coxph(surv_obj ~ alcohol_pw + p25493_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# corticospinal tract
model2_4L <- coxph(surv_obj ~ alcohol_pw + p25496_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_4R <- coxph(surv_obj ~ alcohol_pw + p25497_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# forceps major
model2_5L <- coxph(surv_obj ~ alcohol_pw + p25498_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_5R <- coxph(surv_obj ~ alcohol_pw + p25499_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# inferior fronto-occipital fasciculus
model2_6L <- coxph(surv_obj ~ alcohol_pw + p25500_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_6R <- coxph(surv_obj ~ alcohol_pw + p25501_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# longitudinal fasciculus
model2_7L <- coxph(surv_obj ~ alcohol_pw + p25502_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_7R <- coxph(surv_obj ~ alcohol_pw + p25503_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# medial lemniscus
model2_8L <- coxph(surv_obj ~ alcohol_pw + p25505_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_8R <- coxph(surv_obj ~ alcohol_pw + p25506_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# middle cerebellar peduncle
model2_9 <- coxph(surv_obj ~ alcohol_pw + p25504_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
# parahippocampal part of cingulum
model2_10L <- coxph(surv_obj ~ alcohol_pw + p25494_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_10R <- coxph(surv_obj ~ alcohol_pw + p25495_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# posterior thalamic radiation
model2_11L <- coxph(surv_obj ~ alcohol_pw + p25507_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_11R <- coxph(surv_obj ~ alcohol_pw + p25508_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# superior longitudinal fasciculus
model2_12L <- coxph(surv_obj ~ alcohol_pw + p25509_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_12R <- coxph(surv_obj ~ alcohol_pw + p25510_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# superior thalamic radiation
model2_13L <- coxph(surv_obj ~ alcohol_pw + p25511_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_13R <- coxph(surv_obj ~ alcohol_pw + p25512_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# uncinate fasciculus
model2_14L <- coxph(surv_obj ~ alcohol_pw + p25513_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_14R <- coxph(surv_obj ~ alcohol_pw + p25514_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# 白质纤维束MD指标分析
# acoustic radiation
model2_15L <- coxph(surv_obj ~ alcohol_pw + p25515_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_15R <- coxph(surv_obj ~ alcohol_pw + p25516_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# anterior thalamic radiation
model2_16L <- coxph(surv_obj ~ alcohol_pw + p25517_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_16R <- coxph(surv_obj ~ alcohol_pw + p25518_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# cingulate gyrus part of cingulum
model2_17L <- coxph(surv_obj ~ alcohol_pw + p25519_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_17R <- coxph(surv_obj ~ alcohol_pw + p25520_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# corticospinal tract
model2_18L <- coxph(surv_obj ~ alcohol_pw + p25523_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_18R <- coxph(surv_obj ~ alcohol_pw + p25524_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# forceps major
model2_19L <- coxph(surv_obj ~ alcohol_pw + p25525_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_19R <- coxph(surv_obj ~ alcohol_pw + p25526_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# inferior fronto-occipital fasciculus
model2_20L <- coxph(surv_obj ~ alcohol_pw + p25527_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_20R <- coxph(surv_obj ~ alcohol_pw + p25528_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# longitudinal fasciculus
model2_21L <- coxph(surv_obj ~ alcohol_pw + p25529_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_21R <- coxph(surv_obj ~ alcohol_pw + p25530_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# medial lemniscus
model2_22L <- coxph(surv_obj ~ alcohol_pw + p25532_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_22R <- coxph(surv_obj ~ alcohol_pw + p25533_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# middle cerebellar peduncle
model2_23 <- coxph(surv_obj ~ alcohol_pw + p25531_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri)
# parahippocampal part of cingulum
model2_24L <- coxph(surv_obj ~ alcohol_pw + p25521_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_24R <- coxph(surv_obj ~ alcohol_pw + p25522_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# posterior thalamic radiation
model2_25L <- coxph(surv_obj ~ alcohol_pw + p25534_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_25R <- coxph(surv_obj ~ alcohol_pw + p25535_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# superior longitudinal fasciculus
model2_26L <- coxph(surv_obj ~ alcohol_pw + p25536_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_26R <- coxph(surv_obj ~ alcohol_pw + p25537_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# superior thalamic radiation
model2_27L <- coxph(surv_obj ~ alcohol_pw + p25538_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_27R <- coxph(surv_obj ~ alcohol_pw + p25539_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right
# uncinate fasciculus
model2_28L <- coxph(surv_obj ~ alcohol_pw + p25540_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # left
model2_28R <- coxph(surv_obj ~ alcohol_pw + p25541_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality , data = mri) # right

# 模型数据归纳
summary(model2_1L)
summary(model2_1R)
summary(model2_2L)
summary(model2_2R)
summary(model2_3L)
summary(model2_3R)
summary(model2_4L)
summary(model2_4R)
summary(model2_5L)
summary(model2_5R)
summary(model2_6L)
summary(model2_6R)
summary(model2_7L)
summary(model2_7R)
summary(model2_8L)
summary(model2_8R)
summary(model2_9)
summary(model2_10L)
summary(model2_10R)
summary(model2_11L)
summary(model2_11R)
summary(model2_12L)
summary(model2_12R)
summary(model2_13L)
summary(model2_13R)
summary(model2_14L)
summary(model2_14R)
summary(model2_15L)
summary(model2_15R)
summary(model2_16L)
summary(model2_16R)
summary(model2_17L)
summary(model2_17R)
summary(model2_18L)
summary(model2_18R)
summary(model2_19L)
summary(model2_19R)
summary(model2_20L)
summary(model2_20R)
summary(model2_21L)
summary(model2_21R)
summary(model2_22L)
summary(model2_22R)
summary(model2_23)
summary(model2_24L)
summary(model2_24R)
summary(model2_25L)
summary(model2_25R)
summary(model2_26L)
summary(model2_26R)
summary(model2_27L)
summary(model2_27R)
summary(model2_28L)
summary(model2_28R)


### 模型3 
# 白质纤维束FA指标分析
# acoustic radiation
model3_1L <- coxph(surv_obj ~ alcohol_pw + p25488_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_1R <- coxph(surv_obj ~ alcohol_pw + p25489_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# anterior thalamic radiation
model3_2L <- coxph(surv_obj ~ alcohol_pw + p25490_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_2R <- coxph(surv_obj ~ alcohol_pw + p25491_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# cingulate gyrus part of cingulum
model3_3L <- coxph(surv_obj ~ alcohol_pw + p25492_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_3R <- coxph(surv_obj ~ alcohol_pw + p25493_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# corticospinal tract
model3_4L <- coxph(surv_obj ~ alcohol_pw + p25496_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_4R <- coxph(surv_obj ~ alcohol_pw + p25497_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# forceps major
model3_5L <- coxph(surv_obj ~ alcohol_pw + p25498_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_5R <- coxph(surv_obj ~ alcohol_pw + p25499_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# inferior fronto-occipital fasciculus
model3_6L <- coxph(surv_obj ~ alcohol_pw + p25500_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_6R <- coxph(surv_obj ~ alcohol_pw + p25501_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# longitudinal fasciculus
model3_7L <- coxph(surv_obj ~ alcohol_pw + p25502_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_7R <- coxph(surv_obj ~ alcohol_pw + p25503_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# medial lemniscus
model3_8L <- coxph(surv_obj ~ alcohol_pw + p25505_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_8R <- coxph(surv_obj ~ alcohol_pw + p25506_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# middle cerebellar peduncle
model3_9 <- coxph(surv_obj ~ alcohol_pw + p25504_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
# parahippocampal part of cingulum
model3_10L <- coxph(surv_obj ~ alcohol_pw + p25494_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_10R <- coxph(surv_obj ~ alcohol_pw + p25495_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# posterior thalamic radiation
model3_11L <- coxph(surv_obj ~ alcohol_pw + p25507_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_11R <- coxph(surv_obj ~ alcohol_pw + p25508_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# superior longitudinal fasciculus
model3_12L <- coxph(surv_obj ~ alcohol_pw + p25509_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_12R <- coxph(surv_obj ~ alcohol_pw + p25510_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# superior thalamic radiation
model3_13L <- coxph(surv_obj ~ alcohol_pw + p25511_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_13R <- coxph(surv_obj ~ alcohol_pw + p25512_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# uncinate fasciculus
model3_14L <- coxph(surv_obj ~ alcohol_pw + p25513_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_14R <- coxph(surv_obj ~ alcohol_pw + p25514_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# 白质纤维束MD指标分析
# acoustic radiation
model3_15L <- coxph(surv_obj ~ alcohol_pw + p25515_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_15R <- coxph(surv_obj ~ alcohol_pw + p25516_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# anterior thalamic radiation
model3_16L <- coxph(surv_obj ~ alcohol_pw + p25517_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_16R <- coxph(surv_obj ~ alcohol_pw + p25518_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# cingulate gyrus part of cingulum
model3_17L <- coxph(surv_obj ~ alcohol_pw + p25519_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_17R <- coxph(surv_obj ~ alcohol_pw + p25520_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# corticospinal tract
model3_18L <- coxph(surv_obj ~ alcohol_pw + p25523_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_18R <- coxph(surv_obj ~ alcohol_pw + p25524_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# forceps major
model3_19L <- coxph(surv_obj ~ alcohol_pw + p25525_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_19R <- coxph(surv_obj ~ alcohol_pw + p25526_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# inferior fronto-occipital fasciculus
model3_20L <- coxph(surv_obj ~ alcohol_pw + p25527_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_20R <- coxph(surv_obj ~ alcohol_pw + p25528_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# longitudinal fasciculus
model3_21L <- coxph(surv_obj ~ alcohol_pw + p25529_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_21R <- coxph(surv_obj ~ alcohol_pw + p25530_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# medial lemniscus
model3_22L <- coxph(surv_obj ~ alcohol_pw + p25532_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_22R <- coxph(surv_obj ~ alcohol_pw + p25533_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# middle cerebellar peduncle
model3_23 <- coxph(surv_obj ~ alcohol_pw + p25531_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri)
# parahippocampal part of cingulum
model3_24L <- coxph(surv_obj ~ alcohol_pw + p25521_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_24R <- coxph(surv_obj ~ alcohol_pw + p25522_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# posterior thalamic radiation
model3_25L <- coxph(surv_obj ~ alcohol_pw + p25534_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_25R <- coxph(surv_obj ~ alcohol_pw + p25535_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# superior longitudinal fasciculus
model3_26L <- coxph(surv_obj ~ alcohol_pw + p25536_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_26R <- coxph(surv_obj ~ alcohol_pw + p25537_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# superior thalamic radiation
model3_27L <- coxph(surv_obj ~ alcohol_pw + p25538_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_27R <- coxph(surv_obj ~ alcohol_pw + p25539_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right
# uncinate fasciculus
model3_28L <- coxph(surv_obj ~ alcohol_pw + p25540_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # left
model3_28R <- coxph(surv_obj ~ alcohol_pw + p25541_i2 + preference + with_food + age + sex + ethnicity + education + TDI + smoking + diet_quality + sleep_quality  + blood_pressure + BMI + health_rating + triglycerides + HDL + LDL + glucose + urate, data = mri) # right

# 模型数据归纳
summary(model3_1L)
summary(model3_1R)
summary(model3_2L)
summary(model3_2R)
summary(model3_3L)
summary(model3_3R)
summary(model3_4L)
summary(model3_4R)
summary(model3_5L)
summary(model3_5R)
summary(model3_6L)
summary(model3_6R)
summary(model3_7L)
summary(model3_7R)
summary(model3_8L)
summary(model3_8R)
summary(model3_9)
summary(model3_10L)
summary(model3_10R)
summary(model3_11L)
summary(model3_11R)
summary(model3_12L)
summary(model3_12R)
summary(model3_13L)
summary(model3_13R)
summary(model3_14L)
summary(model3_14R)
summary(model3_15L)
summary(model3_15R)
summary(model3_16L)
summary(model3_16R)
summary(model3_17L)
summary(model3_17R)
summary(model3_18L)
summary(model3_18R)
summary(model3_19L)
summary(model3_19R)
summary(model3_20L)
summary(model3_20R)
summary(model3_21L)
summary(model3_21R)
summary(model3_22L)
summary(model3_22R)
summary(model3_23)
summary(model3_24L)
summary(model3_24R)
summary(model3_25L)
summary(model3_25R)
summary(model3_26L)
summary(model3_26R)
summary(model3_27L)
summary(model3_27R)
summary(model3_28L)
summary(model3_28R)