### Initial processing script for TraceLab TMS Project ###

set.seed(101311)

### Import required packages and functions ###

library(dplyr)
library(ggplot2)
library(brms)
library(tidybayes)
library(emmeans)
library(parameters)
library(modelr)

### Import and preprocess all data ###

if('fulldat.Rdata' %in% dir("./_Scripts/")) {
  source("./_Scripts/0_import.R")
  load("./_Scripts/incomplete_trials.Rdata")
  load("./_Scripts/fulldat.Rdata")
  origins <- points %>%
    group_by(id, session, block, trial) %>%
    summarize(origin.x = x[1], origin.y = y[1])
} else {
  source("./_Scripts/1_preprocessing.R")
}

#create visualization directory if needed
if('_Vis' %in% dir("./")){
}else{
  dir.create("./_Vis")
}

### Get initial sample info ###

init_id_info <- fulldat %>%
  group_by(participant) %>%
  summarize(
    db_id = db_id[1],
    sex = sex[1],
    age = age[1],
    handedness = handedness[1],
    trialcount = n()
  )

### Filter out problem participants ###

# Drop IDs who did not complete enough useable physical trials

phys_trial_counts <- fulldat %>%
  group_by(db_id, session_num) %>%
  summarize(phys_count = sum(!is.na(mt_clip)))

incomplete_ids <- subset(phys_trial_counts, phys_count < 16)$db_id
a <- subset(fulldat, !(db_id %in% incomplete_ids))

### Get final sample info ###

id_info <- subset(init_id_info, db_id %in% unique(a$db_id))
group_info <- group_by(id_info) %>%
  summarise(
    age_mean = mean(age),
    age_sd = sd(age),
    n_female = sum(sex=='f'),
    n_left = sum(handedness=='l')
  )

### Additional data cleaning ###

# Drop all physical trials with unusable tracing data

a <- subset(a, !(is.na(mt_clip)))

#look at all the incomplete trials

merge_key <- c(
  "db_id" = "id",
  "session_num" = "session",
  "block_num" = "block",
  "trial_num" = "trial"
)

a_inc <- right_join(fulldat, incomplete_trials, by = merge_key) %>%
  select(
    c(db_id, session_num, block_num, trial_num, figure_type, sinuosity,
      totabscurv, turnangle_sum, turnangle_mean, turnangle_sd,
      approx_en, sample_en
    )
  )

# Drop unnecessary columns - update once columns clear up

a2 <- a %>%
  select(-c(participant, age, feedback_type, figure_set)) %>%
  select(-c(control_question, control_response, figure_file, trace_file)) %>%
  rename(
    id = db_id, session = session_num,
    block = block_num, trial = trial_num
  ) %>%
  select(id, session, block, trial, sex, handedness, everything())

# Fix incorrect stimulus_mt values

true_stim_mt <- frames %>%
  group_by(id, session, block, trial) %>%
  summarize(actual_mt = max(time))

a2 <- a2 %>%
  left_join(true_stim_mt, by = c("id", "session", "block", "trial")) %>%
  mutate(stimulus_mt = actual_mt) %>%
  select(-actual_mt)

### Do some descriptives ###

a3 <- subset(a2, abs(stimulus_mt - mt_clip) < 2)

# Standardize / transform / factorize variables before modelling

# pick outcome variable
outcome_var='dtw_proc_err_mean'
#eqd_err_mean.x = PROC, eqd_err_mean.y = EQD, eqd_err_mean.x.x = EQD&PROC, dtw_err_mean = DTW, eqd_err_mean.y.y. = DTW&PROC

#remove equidistant function that fails
if(outcome_var == c('dtw_proc_err_mean')){
  a4 <- filter(a3,!is.na(dtw_proc_err_mean))
}
#transform
a4 <- a4 %>%
  mutate(count=0:(n()-1), .before = 5) %>%
  group_by(id, session,learning_type) %>%
  mutate(exposure=0:(n()-1),.before = 6) %>%
  group_by(id, session)

a4$log_dtw_mean <- log(a4$dtw_proc_err_mean)
a4$log_stim_vel <- log(a4$real_length/(a4$stimulus_gt/1000))

#standardize
a4$log_dtw_mean_z <- (a4$log_dtw_mean - mean(a4$log_dtw_mean)) / sd(a4$log_dtw_mean)
a4$log_stim_vel_z <- (a4$log_stim_vel-mean(a4$log_stim_vel))-sd(a4$log_stim_vel)
a4$turnangle_sum_z <- (a4$turnangle_sum - median(a4$turnangle_sum)) / sd(subset(a4, figure_type == "random")$turnangle_sum)

#create factors
a4$id <- as.factor(a4$id)
a4$stim_duration <- factor(a4$stimulus_gt, levels = c("500", "1000", "1500", "2000", "2500"), ordered = TRUE)
a4$learning_type <- as.factor(a4$learning_type)

#set contrast coding
contrasts(a4$learning_type) <- matrix(c(-0.5,0.5,0,-1/3,-1/3,2/3),ncol = 2)

#plot descriptives

pd <- position_dodge(0.3)
makinStuffPretty <-   theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    #legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 22),
    axis.text = element_text(colour = "black",size = 22),
    axis.ticks = element_line(colour = "black",size = 1),
  )

alldat <- ggplot(filter(a4,exposure<29),aes(y=dtw_proc_err_mean,x=exposure+1,colour=learning_type))+
  #eqd_err_mean.x = PROC, eqd_err_mean.y = EQD, eqd_err_mean.x.x = EQD&PROC, dtw_err_mean = DTW, eqd_err_mean.y.y = DTW&PROC
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = position_dodge(0.5)) +
  makinStuffPretty +
  scale_colour_manual(values = c("grey80","grey40","black")) +
  labs(colour = "Learning Type", x="Exposure",y="Mean Trajectory Error (px)")

show(alldat)

ggsave(filename="alldat.png",
       path="./_Vis/",
       plot=alldat,
       dpi=600,
       width=30,
       height=15,
       units="cm")

#plot first 5 trials

first5 <- ggplot(subset(a4, exposure <5),aes(y=dtw_proc_err_mean,x=exposure+1,colour=learning_type))+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange", position = pd) +
  makinStuffPretty +
  scale_colour_manual(values = c("grey80","grey40","black")) +
  labs(colour = "Learning Type", x="Exposure",y="Mean Trajectory Error (px)")

show(first5)

ggsave(filename="first5.png",
       path="./_Vis/",
       plot=first5,
       dpi=600,
       width=20,
       height=15,
       units="cm")

### Modeling ###

model <- log_dtw_mean_z ~  exposure * learning_type * log_stim_vel_z +
  turnangle_sum_z +
  (1+log_stim_vel_z * exposure * learning_type | id)

altmod <- brm(model,
                data = filter(a4,exposure<30),
                prior = c(
                  prior(normal(0, 2), class = b),
                  prior(exponential(1), class = sd),
                  prior(exponential(1), class = sigma)
                ),
                iter = 20000, warmup = 5000, chains = 8, cores = 8, control = list(adapt_delta=.99, max_treedepth = 20)
)

res <- tibble(parameters(altmod, ci = 0.89, ci_method = "hdi",test=c("pd","rope"),rope_ci = 0.95))

# Generate draws from model for visualizing linear predictors ###

# Generate draws across values of stimulus animation velocity, holding complexity at
# repeated mean and plot

altmod_fitted <- ungroup(a4) %>%
  data_grid(
    turnangle_sum_z = 0,
    exposure=median(a4$exposure),
    learning_type=unique(a4$learning_type),
    log_stim_vel_z = seq_range(log_stim_vel_z, n = 101)
  ) %>%
  add_linpred_draws(altmod, re_formula = NA, ndraws = 2000) %>%
  mutate(
    log_err = (.linpred * sd(a4$log_dtw_mean)) + mean(a4$log_dtw_mean),
    shape_err = exp(log_err),
    log_stim_vel = (log_stim_vel_z * sd(a4$log_stim_vel)) + mean(a4$log_stim_vel),
    stim_vel = exp(log_stim_vel)
  )

speedXlearning <- ggplot(altmod_fitted, aes(x = stim_vel, y = shape_err,group = learning_type, color = learning_type)) +
  stat_lineribbon(.width = .90, point_interval = median_qi,show.legend = TRUE,
                  fill=alpha(c("grey"), 0.25)) +
  xlab("Stimulus Animation Velocity (px/s)") +
  ylab("Mean Trajectory Error (px)") +
  scale_colour_manual(values = c("grey80","grey40","black")) +
  makinStuffPretty +
  guides(color = guide_legend(order = 1),
         fill = "none")

show(speedXlearning)


ggsave(filename="speedXlearning.png",
       path="./_Vis/",
       plot=speedXlearning,
       dpi=600,
       width=20,
       height=15,
       units="cm")

# Generate draws across values of turnangle_sum_z2

altmod_fitted <- ungroup(a4) %>%
  data_grid(
    turnangle_sum_z = seq_range(turnangle_sum_z, n = 5),
    exposure=median(a4$exposure),
    learning_type=unique(a4$learning_type),
    log_stim_vel_z = median(log_stim_vel_z)
  ) %>%
  add_linpred_draws(altmod, re_formula = NA, ndraws = 2000) %>%
  mutate(
    log_err = (.linpred * sd(a4$log_dtw_mean)) + mean(a4$log_dtw_mean),
    shape_err = exp(log_err),
    log_stim_vel = (log_stim_vel_z * sd(a4$log_stim_vel)) + mean(a4$log_stim_vel),
    stim_vel = exp(log_stim_vel),
    turnangle = (turnangle_sum_z * sd(subset(a4, figure_type == "random")$turnangle_sum)) + median(a4$turnangle_sum),
  )

compAll <- ggplot(altmod_fitted, aes(x = turnangle, y = shape_err)) +
  stat_lineribbon(colour="black",.width = .90, point_interval = median_qi,
                  fill=alpha(c("grey"), 0.25)) +
  xlab("Total Absolute Curvature (rad)") +
  ylab("Mean Trajectory Error (px)") +
  makinStuffPretty

show(compAll)


ggsave(filename="compAll.png",
       path="./_Vis/",
       plot=compAll,
       dpi=600,
       width=20,
       height=15,
       units="cm")

# Restrict Analysis to Block 1 #

a5 <- filter(a4)%>%
  group_by(id,learning_type) %>%
  filter(id!=41) %>%
  ungroup()

contrasts(a5$learning_type) <- matrix(c(-0.5,0.5,0,-1/3,-1/3,2/3),ncol = 2)

altmod1 <- brm(model,
                data = filter(a5,exposure <5),
                prior = c(
                  prior(normal(0, 2), class = b),
                  prior(exponential(1), class = sd),
                  prior(exponential(1), class = sigma)
                ),
                iter = 20000, warmup = 5000, chains = 8, cores = 8,
                control = list(max_treedepth = 20, adapt_delta=.99)
)

res5 <- tibble(parameters(altmod1, ci = 0.89, ci_method = "hdi",test=c("pd","rope"),rope_ci = 0.95))

# interaction of speed*exposure*figure
altmod1_fitted <- ungroup(a5) %>%
  data_grid(
    turnangle_sum_z = 0,
    exposure=c(0:4),
    learning_type=unique(a5$learning_type),
    log_stim_vel_z = seq_range(log_stim_vel_z, n = 101)
  ) %>%
  add_linpred_draws(altmod1, re_formula = NA, ndraws = 2000) %>%
  mutate(
    log_err = (.linpred * sd(a5$log_dtw_mean)) + mean(a5$log_dtw_mean),
    shape_err = exp(log_err),
    log_stim_vel = (log_stim_vel_z * sd(a5$log_stim_vel)) + mean(a5$log_stim_vel),
    stim_vel = exp(log_stim_vel),
    exposure=exposure+1
  )

speedXexpXlearning <- ggplot(altmod1_fitted, aes(x = stim_vel, y = shape_err,group = learning_type, color = learning_type)) +
  stat_lineribbon(.width = .90, point_interval = median_qi,show.legend = TRUE,
                  fill=alpha(c("grey"), 0.25)) +
  labs(group = "Trial Type", color = "Trial Type") +
  xlab("Stimulus Animation Velocity (px/s)") +
  ylab("Mean Trajectory Error (px)") +
  facet_wrap(~exposure, 
             nrow=1,
             labeller = as_labeller(c("1" = "Exposure 1", "2" = "Exposure 2", "3" = "Exposure 3","4" = "Exposure 4", "5" = "Exposure 5")))+
  scale_colour_manual(values = c("grey80","grey40","black")) +
  makinStuffPretty +
  theme(
    axis.text.x = element_text(angle=90,vjust=.5)
  ) +
  guides(color = guide_legend(order = 1),
         fill = "none")

show(speedXexpXlearning)

ggsave(filename="speedXexpXlearning.png",
       path="./_Vis/",
       plot=speedXexpXlearning,
       dpi=600,
       width=30,
       height=15,
       units="cm")

# main effect of complexity
altmod1_fitted <- ungroup(a5) %>%
  data_grid(
    turnangle_sum_z = seq_range(turnangle_sum_z, n = 101),
    exposure=c(1:4),
    learning_type=unique(a5$learning_type),
    log_stim_vel_z = median(log_stim_vel_z)
  ) %>%
  add_linpred_draws(altmod1, re_formula = NA, ndraws = 2000) %>%
  mutate(
    log_err = (.linpred * sd(a5$log_dtw_mean)) + mean(a5$log_dtw_mean),
    shape_err = exp(log_err),
    log_stim_vel = (log_stim_vel_z * sd(a5$log_stim_vel)) + mean(a5$log_stim_vel),
    stim_vel = exp(log_stim_vel),
    turnangle = (turnangle_sum_z * sd(subset(a4, figure_type == "random")$turnangle_sum)) + median(a4$turnangle_sum)
  )

comp5 <- ggplot(altmod1_fitted, aes(x = turnangle, y = shape_err)) +
  stat_lineribbon(colour="black",.width = .90, point_interval = median_qi,
                  fill=alpha(c("grey"), 0.25)) +
  xlab("Total Absolute Curvature (rad)") +
  ylab("Mean Trajectory Error (px)") +
  makinStuffPretty

show(comp5)

ggsave(filename="comp5.png",
       path="./_Vis/",
       plot=comp5,
       dpi=600,
       width=20,
       height=15,
       units="cm")
