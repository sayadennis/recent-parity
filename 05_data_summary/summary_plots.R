library(tidyverse)

set.seed(123)

ds <- data.frame(
  ID = seq.int(60),
  situation = rep(c("Soccer", "Basketball", "snooker"), 20),
  attitude = c("aggressive", "assertive", "neutral"),
  response_t1 = c("wrong", "Right", "I Don't"),
  response_t2 = rep(c("wrong", "Right", "I Don't"), times = c(10, 35, 15))
)

ds %>%
  gather(key = "Time", value, response_t1:response_t2) -> j
#> Warning: attributes are not identical across measure variables;
#> they will be dropped

j %>%
  group_by(attitude, Time, situation, value) %>%
  summarise(n = n()) %>%
  ggplot(., aes(x = value, y = n, fill = Time)) +
  geom_col(width = 0.5, position = position_dodge(preserve = "single", width = 0.6)) +
  facet_wrap(~situation * attitude, ncol = 3)

