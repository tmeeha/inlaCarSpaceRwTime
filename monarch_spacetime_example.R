# set up -----------------------------------------------------------------------
setwd("C:/Users/tmeehan/Box/-.Tim.Meehan individual/monarch_space_time")

# load some libraries
library(ggplot2); library(cowplot); library(bbsBayes); library(tidybayes)
library(ggpattern); library(coda); library(spdep); library(brinla); library(scico)
library(inlatools); library(sf); library(INLA); library(tidyr); 
library(stringr); library(dplyr)

# set some options
inla.setOption(inla.mode="experimental")
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors=F)

# helper for getting model scores
get_bpic <- function(mod=res1, pred_idx=pred_idx1){
  require(INLA)
  cpo_vec <- mod$gcpo$gcpo[-pred_idx]
  nas <- which(is.na(cpo_vec))
  if(any(is.na(cpo_vec))) 
    warning(paste("Some CPO values are NA. BPIC may not be", 
                  "reliable for model comparison.")) 
  out1 <- as.numeric(-2*sum(log(cpo_vec), na.rm=T))
  out2 <- as.numeric(mean(log(cpo_vec), na.rm=T))
  out3 <- as.integer(nas)
  out <- list(`-2*sum(log(cpo_vec), na.rm=T)`=out1,
              `mean(log(cpo_vec), na.rm=T)`=out2,
              `number_na`=length(out3),
              `na_indices`=out3)
  return(out)
}

# helper for getting model predictions
get_prediction_samples <- function(mod_obj=res1, post_samp_size=1000,
                                   pred_idx=pred_idx1, pred_data=d3,
                                   pred_id_cols=c("strat_name", "year"),
                                   out_form=c("tidy_tibble", "mcmc_list")
                                   ){
  require(coda); require(tidyr); library(tidybayes); require(INLA); 
  require(dplyr)
  s1 <- INLA::inla.posterior.sample(post_samp_size, mod_obj)
  s2 <- INLA::inla.posterior.sample.eval("Predictor", s1)[pred_idx,] 
  s2 <- t(exp(s2))
  ids1 <- pred_data %>% dplyr::slice(pred_idx) %>% 
    dplyr::select(pred_id_cols)  
  colnames(s2) <- paste0("index[", paste(ids1$strat_name, 
                                                ids1$year, sep=","), "]")
  s3 <- as.mcmc.list(mcmc(s2))
  if(out_form=="tidy_tibble"){
    s3 <- spread_draws(s3, index[strat_name, year])
  }
  return(s3)
}
# ------------------------------------------------------------------------------


# maps and data ----------------------------------------------------------------
# bring in bbs stratum map from bbsBayes package
map1 <- load_map(stratify_by="bbs_usgs") %>% rename_all(tolower) %>% 
  select(2, 3) %>% rename(strat_name=st_12, strat_A=area_1)

# choose the year and date range using some exploratory plots
start_yr <- 1994
read.csv("monarch_site_summary.csv") %>% 
  filter(between(j_date, 152, 243),
         year>=start_yr) %>% 
  distinct(year, site) %>% group_by(year) %>% 
  summarise(n=length(unique(site))) %>% 
  ggplot() + geom_line(aes(x=year, y=n))
read.csv("monarch_site_summary.csv") %>% 
  filter(between(j_date, 152, 243), # jdate relates to site latitude
         year>=start_yr) %>% # so don't statistically control for it
  mutate(lat=as.numeric(str_split(site, "_", simplify=T)[,3])) %>% 
  group_by(j_date) %>% 
  summarise(n=mean(n_monarch),
            l=mean(lat)) %>% 
  ggplot() + geom_point(aes(x=j_date, y=n, col=l, size=l)) +
  scale_color_distiller(palette="Spectral")
read.csv("monarch_site_summary.csv") %>% 
  filter(between(j_date, 152, 243), # timing has been consistent over years
         year>=start_yr) %>% 
  group_by(year) %>% 
  summarise(n=mean(j_date)) %>% 
  ggplot() + geom_line(aes(x=year, y=n))
read.csv("monarch_site_summary.csv") %>% 
  filter(between(j_date, 152, 243), # hours have increased so statistically
         year>=start_yr) %>% # adjust for them
  group_by(year) %>% 
  summarise(n=mean(hours)) %>% 
  ggplot() + geom_line(aes(x=year, y=n))

# load naba monarch data. this data is ALREADY ZERO FILLED.
# if you data is not already zero-filled, you should do it now.
# what does this mean? many count databases are set up to only
# store records when a critter is observed and ignore zeros for
# efficiency. but if a count was conducted and zero critters 
# were observed, then you want this information in a trend 
# analysis. you need to know when all counts occurred, and the 
# effort involved, and then you can fill in zeros where a count 
# was conducted but no critters were seen. this is an annoying
# but necessary data munging step.
d0 <- read.csv("monarch_site_summary.csv") %>% 
  filter(between(j_date, 180, 220),
         year>=start_yr) %>%  
  select(site, year, n_hours=hours, n_monarch) %>% 
  group_by(site) %>% 
  mutate(n_detects=sum(n_monarch>0),
         n_years=length(unique(year)),
         first_yr=min(year),
         last_yr=max(year),
         yr_span=last_yr-first_yr+1,
         detect_prop=n_detects/yr_span,
         lon=str_split(site, "_", simplify = T)[,2],
         lat=str_split(site, "_", simplify = T)[,3]) %>% 
  ungroup()

# create a map of all sites, regardless of monarch presence.
# this will be used to calculate z, the proportion of sites 
# in a stratum with monarchs detected
d0 %>% filter(year %in% seq(min(year), max(year), 2)) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=st_crs(map1)) %>% 
  ggplot() + geom_sf(aes(col=sqrt(n_monarch/n_hours)), size=1) + 
  facet_wrap(~year) + scale_color_distiller(palette="Spectral")
sites1 <- d0 %>% distinct(site, lon, lat) %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=st_crs(map1))

# pair each site with its respective bbs stratum
within1 <- sapply(st_within(sites1, map1), function(z) if (length(z)==0) 
  NA_integer_ else z[1])
near1 <- st_nearest_feature(sites1, map1) 
poly_per_site <- ifelse(is.na(within1), near1, within1)
site_strat_key <- sites1 %>% 
  mutate(strat_name=unique(map1$strat_name)[poly_per_site]) %>% 
  st_drop_geometry() %>% 
  right_join(map1 %>% st_drop_geometry()) %>% 
  select(site, strat_name, strat_A)

# join with d0 to get the number of naba sites and z per stratum
strat_table <- full_join(d0, site_strat_key) %>% 
  group_by(strat_name, site) %>%
  summarise(monarch_detected=as.numeric(sum(n_monarch, na.rm=T)>0),
            strat_A=mean(strat_A, na.rm=T)) %>% 
  ungroup() %>% 
  group_by(strat_name) %>% 
  summarise(n_sites=length(unique(site[!is.na(site)])),
            strat_A=mean(strat_A),
            strat_z=sum(monarch_detected)/n_sites) %>% 
  ungroup() %>% mutate(strat_z=ifelse(is.na(strat_z), 0, strat_z))

# select sites for modeling based on the number of detections
# and merge strata info, z, and A.
# choosing 5 means that a counts can be included if there are only 5 years
# of counts at a site, but a monarch was detected each year, or if a count
# site recorded monarchs on 5 out of 23 years
names(d0)
d1 <- d0 %>% left_join(site_strat_key) %>% left_join(strat_table) %>% 
  filter(n_detects>=5) %>% 
  group_by(strat_name) %>% 
  mutate(n_sites_modeled=length(unique(site[!is.na(site)]))) %>% 
  ungroup() %>% 
  filter(n_sites_modeled>=3) %>% 
  select(-n_detects, -n_years, -first_yr, -last_yr, -yr_span, -detect_prop)
unique(d1$site)
summary(d1)

# create a map of the sites selected for modeling
d1 %>% filter(year %in% seq(min(year), max(year), 2)) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=st_crs(map1)) %>% 
  ggplot() + geom_sf(aes(col=sqrt(n_monarch/n_hours)), size=1) + 
  facet_wrap(~year) + scale_color_distiller(palette="Spectral")
sites2 <- d1 %>% group_by(site, lon, lat) %>% 
  summarise(n_years=sum(year>1900)) %>%
  ungroup() %>%
  mutate(n_years=ifelse(n_years > length(unique(d1$year)), 
                        length(unique(d1$year)), n_years)) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=st_crs(map1))

# slim down the bbs stratum map to strata with selected sites for modeling.
# this includes all contiguous strata regardless of count site.
# keeping the empty strata retains spatial relationships between 
# polygons. dropping them distorts space.
# this creates map2 which has indices that will be used in neighborhood
# generation, below.
area1 <- st_convex_hull(st_combine(sites2)) %>% st_as_sf()
strat_in_area <- sapply(st_intersects(map1, area1), 
                        function(z) if (length(z)==0) 
  NA_integer_ else z[1])
map2 <- map1[which(strat_in_area==1), ] %>% 
  mutate(strat_idx1=as.integer(factor(strat_name))) %>% 
  select(strat_name, strat_idx1, strat_A) %>% 
  left_join(strat_table) %>% 
  arrange(strat_name)
ggplot() + geom_sf(data=map2) + geom_sf(data=sites2, aes(size=n_years)) + 
  geom_sf(data=area1, col="darkblue", fill=NA) + theme_bw()

# pair the selected sites with the included stratum IDs.
# this pairing is done to give sites the new stratum IDs from map2.
within2 <- sapply(st_within(sites2, map2), function(z) if (length(z)==0)
  NA_integer_ else z[1])
near2 <- st_nearest_feature(sites2, map2)
poly_per_site <- ifelse(is.na(within2), near2, within2)
site_strat_key <- sites2 %>%
  mutate(strat_name=unique(map2$strat_name)[poly_per_site]) %>%
  select(site, strat_name, everything()) %>% st_drop_geometry() %>%
  right_join(map2 %>% st_drop_geometry())

# join with d1 and filter out strata with < 3 sites.
# this is done so that a stratum trend is not based on 1 or 2 sites.
# this means a trend will be generated from as little as 3 sites (below) 
# times 5 detections (above).
d2 <- d1 %>% full_join(site_strat_key) %>% 
  select(strat_name, strat_idx1, strat_A, strat_z, year, everything()) %>% 
  mutate(log_hrs=log(n_hours),
         site_idx1=as.integer(factor(site)),
         year_idx1=as.integer(factor(year))) %>% 
  arrange(strat_idx1, year_idx1, site_idx1) %>% 
  filter(!is.na(year))
d3 <- d2 %>% select(strat_idx1, year_idx1, 
                    site_idx1, log_hrs, count=n_monarch)

# create prediction rows per stratum and year. these are rows with mean hours
# effort and ignoring site intercepts (NA), and setting count to NA so that 
# INLA will give it a prediction.
year_strat_combs <- d3 %>% tidyr::expand(year_idx1, strat_idx1)
pred_rows <- year_strat_combs %>% 
  left_join(d3) %>% 
  mutate(log_hrs=log(mean(exp(d2$log_hrs))),
         site_idx1=NA, count=NA) %>% 
  distinct(strat_idx1, year_idx1, .keep_all=T)

# add prediction rows to dataset
d3 <- d3 %>% bind_rows(pred_rows) %>% mutate(strat_idx2=strat_idx1) %>%
  st_drop_geometry() %>% 
  left_join(map2 %>% select(strat_name, strat_idx1)) %>% 
  st_drop_geometry() %>% 
  mutate(year=year_idx1+min(d2$year)-1) %>% 
  select(strat_name, year, log_hrs, count, site_idx1, strat_idx2, 
         strat_idx1, year_idx1)
summary(d3)
pred_idx1 <- (nrow(d2)+1):nrow(d3)

# make the neighborhood matrix and graph for INLA using map2
ggplot() + geom_sf(data=map2) + 
  geom_sf(data=sites2 %>% filter(site %in% d2$site), aes(size=n_years, col=n_years)) + 
  scale_color_scico(palette = "vik", direction=-1,
                     begin=0.4, end=0.1) +
  geom_sf(data=area1, col="black", lty=2, fill=NA) + theme_bw()
nb1 <- poly2nb(as(map2, "Spatial"))
summary(nb1)
nb2INLA("map.adj", nb1)
g1 <- inla.read.graph(filename = "map.adj")
plot(g1)
# ------------------------------------------------------------------------------


# model and output -------------------------------------------------------------
# set the priors for the model components
bym_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                  phi = list(prior = "pc", param = c(0.5, 0.5)))
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

# make the space-time model formula
form <- count ~ 1 + log_hrs +
  f(site_idx1, model="iid", constr=T, hyper=prec_prior) + # rand site intercept
  f(strat_idx1, log_hrs, model="besag", graph=g1, constr=T, scale.model=T, 
    hyper=prec_prior) + # svc effort correction
  f(strat_idx2, model="bym2", graph=g1, constr=T, #space-time part
    scale.model=T, hyper=bym_prior, 
    group=year_idx1, control.group=list(model="rw1", hyper=prec_prior)) 

# # what about a more-flexible effort correction?
# form <- count ~ 1 + 
#   f(site_idx1, model="iid", constr=T, hyper=prec_prior) + # rand site intercept
#   f(round(log_hrs, 1), model="rw2", constr=T, scale.model=T, hyper=prec_prior,
#     group=strat_idx1, control.group=list(model="iid", hyper=prec_prior)) + # svc effort correction
#   f(strat_idx2, model="bym2", graph=g1, constr=T, #space-time part
#     scale.model=T, hyper=bym_prior, 
#     group=year_idx1, control.group=list(model="rw1", hyper=prec_prior)) 

# run the space-time model
res1 <- inla(form, data=d3,
            family = "nbinomial", 
            control.predictor=list(compute=T, link=1),
            control.compute=list(config=T, waic=T, cpo=T,
                                 control.gcpo=list(enable=T)),
            verbose=T, inla.mode="experimental",
            control.inla = list(strategy = 'adaptive', int.strategy = 'eb'))

# view some model fit diagnostics
fp1 <- inlatools::dispersion_check(res1) %>% plot() + 
  theme_bw() + theme(plot.title=element_blank()) +
  labs(x="Dispersion", y="Density") 
fp2 <- inlatools::fast_distribution_check(res1) %>% 
  plot() + 
  theme_bw() + 
  labs(x="Predicted abundance index", y="Observed / predicted")
fp3 <- ggplot(data.frame(x=res1$cpo$pit[-pred_idx1])) +
  geom_histogram(aes(x=x), col="white", fill=scico(10, palette = 'vik')[4]) + 
  theme_bw() +
  labs(x="PIT value", y="Count")
fp4 <- ggplot(data.frame(predicted=res1$summary.fitted.values$mean, 
       observed=d3$count), aes(x=sqrt(observed), y=sqrt(predicted))) + 
  geom_hex() + labs(x="Square root observed", y="Square root predicted") +
  scale_fill_scico("Number of\nobservations", palette = "vik", midpoint=200, 
                   begin=0.1, end=0.9) +
  theme_bw() + 
  geom_abline(aes(intercept=0, slope=1), col="black", lty=2)
plot_grid(fp1, fp2, fp3, fp4)
cor(res1$summary.fitted.values$mean, d3$count, use="complete.obs")

# view hyperparameter and fixed effect summary
bri.hyperpar.summary(res1)[-1,c(1,3,5)]
res1$summary.fixed[,c(1,3,5)]
summary(res1$summary.random$strat_idx1$mean)
summary(res1$summary.random$strat_idx2$mean)

# get some model scores in case you want to do model selection
bpic1 <- get_bpic(res1, pred_idx1); bpic1$`-2*sum(log(cpo_vec), na.rm=T)`
waic1 <- res1$waic$waic; waic1

# extract predicted abundance indices from the prediction rows
# added above to d3. then multiply by stratum scalars A and z
d3$st_fit <- exp(res1$summary.linear.predictor$mean)
d3$st_lcl <- exp(res1$summary.linear.predictor$`0.025quant`)
d3$st_ucl <- exp(res1$summary.linear.predictor$`0.975quant`)
pred1 <- d3 %>% slice(pred_idx1) %>% 
  left_join(map2) %>% 
  left_join(d2 %>% select(strat_name, n_sites_modeled) %>% 
              distinct()) %>% 
  st_as_sf() %>% 
  select(strat_name, strat_A, strat_z, n_sites, n_sites_modeled, 
         year, strat_idx1, year_idx1, st_fit, st_lcl, st_ucl,
         geometry) %>% ungroup() %>% 
  mutate(st_fit=st_fit*strat_A*strat_z,
         st_lcl=st_lcl*strat_A*strat_z,
         st_ucl=st_ucl*strat_A*strat_z)

# map abundance indices per stratum per 1 year
ggplot() + 
  geom_sf(data=pred1 %>% 
            filter(year %in% seq(min(year), max(year), by=1)), 
          aes(fill = st_fit/1000000)) +
  facet_wrap(~year, dir = "h", ncol = 7) +
  ggtitle("Posterior mean BBS abundance index") + theme_bw() +
  theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) +
  scale_fill_scico("Relative\nabundance", palette = "vik", midpoint=4, 
                   begin=0.1, end=0.9)

# plot a selection of stratum trajectories 
ggplot(pred1 %>% filter(grepl('ON|NJ|NY|PA', strat_name)), 
       aes(x=year, y=st_fit/1000000, ymin=st_lcl/1000000, 
           ymax=st_ucl/1000000)) +
  geom_ribbon(alpha=0.4) +
  geom_line(col="darkblue") +
  facet_wrap(~strat_name, scales="free") +
  labs(x="Year", y="Spacetime fit (95% CrI)") +
  theme_bw() + scale_x_continuous(breaks=seq(1994, 2018, 4))
  
# # what about a more-flexible effort correction?
# res1$summary.random$`round(log_hrs, 1)`[,1:2] %>% 
#   mutate(idx=rep(1:169, each=43)) %>% 
#   filter(idx==sample(1:169, 20)) %>% 
#   ggplot() + geom_line(aes(x=ID, y=mean)) +
#   facet_wrap(~idx, scales="free")
  
# get a posterior sample of abundance indices. sample size is initially 
# set low for model development. it should be around 5000 for the final run.
samps1 <- get_prediction_samples(mod_obj=res1, 
                                   post_samp_size=1000,
                                   pred_idx=pred_idx1,
                                   pred_data=d3,
                                   pred_id_cols=c("strat_name", "year"),
                                   out_form="tidy_tibble")
samps1 <- samps1 %>% filter(strat_name %in% pred1$strat_name)
sum(samps1$index)

# # weight by land use change
# samps1 <- samps1 %>% 
#   mutate(bcr=str_split(strat_name, "-", simplify=T)[,3]) %>% 
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year<2000, index/0.5, index)) %>% 
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2001, index/0.55, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2002, index/0.6, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2003, index/0.65, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2004, index/0.7, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2005, index/0.75, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2006, index/0.8, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2007, index/0.85, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year==2008, index/0.9, index)) %>%
#   mutate(index=ifelse((bcr=="22" | bcr=="22") & year>2008, index/0.95, index)) %>%
#   select(-bcr)
# sum(samps1$index)
# ------------------------------------------------------------------------------


# aggregation and trends bcr ---------------------------------------------------
# extract the areas and z proportions per stratum so scaling 
# can be done to samples
Azs <-  pred1 %>% st_drop_geometry() %>% 
  distinct(strat_name, strat_A, strat_z) %>% 
  mutate(strat_Az=strat_A*strat_z)
samps3 <- Azs %>% 
  right_join(samps1) %>% 
  mutate(country=str_split(strat_name, "-", simplify=T)[,1],
         state=str_split(strat_name, "-", simplify=T)[,2],
         bcr=str_split(strat_name, "-", simplify=T)[,3],
         strat_nAz=strat_Az*index) %>% 
  group_by(bcr, year, .draw) %>% 
  summarise(aggregate_N=sum(strat_nAz), aggregate_A=sum(strat_A)) %>% 
  ungroup()

# trajectories per bcr
samps3 %>% 
  ggplot(aes(x=year, y=aggregate_N/1000000)) +
      stat_lineribbon(.width=c(0.5, 0.8, 0.95), lwd=0.8) + 
  scale_fill_scico_d(palette = "vik", 
                     begin=0.45, end=0.15) +
      facet_wrap(~bcr, scale="free") + theme_bw()

# poisson regression trends per bcr since 2008
cbc_bcr_reg_trends <- samps3 %>% 
  filter(year>=2008) %>% 
  group_by(bcr, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(bcr) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))

# map for mapping regression trends by bcr
map3 <- load_map(stratify_by="bcr") %>% rename_all(tolower) %>% 
  rename(bcr=st_12, strat_A=area_sq_km) %>%
  mutate(bcr=gsub("BCR","", bcr)) %>% 
  filter(bcr %in% cbc_bcr_reg_trends$bcr) %>% 
  right_join(cbc_bcr_reg_trends) 

# map of the trends by bcr since 2008
p2008 <- ggplot() + 
  geom_sf(data=map3, aes(fill=mean_trend)) + 
  geom_sf_pattern(data=map3 %>% filter(diff_zero==1),
                  fill=NA,
                  pattern='stripe',
                  pattern_density=0.5,
                  pattern_fill=NA) +
  geom_sf(data=sites2, alpha=0.4, aes(size=n_years), pch=16) +
  scale_size_continuous("number\nyears", range=c(0.2, 2), guide="none") +
  scale_fill_scico("Trend\n2008-\n2018", palette = "vik", midpoint=0, 
                   begin=0.1, end=0.9) + 
  theme_bw(); p2008

# poisson regression trends per bcr for all years
cbc_bcr_reg_trends <- samps3 %>% 
  group_by(bcr, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(bcr) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))

# map regression trends
map4 <- load_map(stratify_by="bcr") %>% rename_all(tolower) %>% 
  rename(bcr=st_12, strat_A=area_sq_km) %>%
  mutate(bcr=gsub("BCR","", bcr)) %>% 
  filter(bcr %in% cbc_bcr_reg_trends$bcr) %>% 
  right_join(cbc_bcr_reg_trends) 

# map it all together
pall <- ggplot() + 
  geom_sf(data=map4, aes(fill=mean_trend)) + 
  geom_sf_pattern(data=map4 %>% filter(diff_zero==1),
                  fill=NA,
                  pattern='stripe',
                  pattern_density=0.5,
                  pattern_fill=NA) +
  geom_sf(data=sites2, alpha=0.5, aes(size=n_years), pch=16) +
  scale_size_continuous("number\nyears", range=c(0.2, 2), guide="none") +
  scale_fill_scico(paste0("Trend\n",start_yr,"-\n2018"), palette = "vik", 
                   midpoint=0, begin=0.1, end=0.9) + 
  theme_bw(); pall
plot_grid(pall, p2008)
# ------------------------------------------------------------------------------


# aggregations and trends country ----------------------------------------------
# aggregate the abundance indices across stratum and sample per country
# # by midwest
# samps4 <- Azs %>%
#   right_join(samps1) %>%
#   mutate(country=str_split(strat_name, "-", simplify=T)[,1],
#          state=str_split(strat_name, "-", simplify=T)[,2],
#          bcr=str_split(strat_name, "-", simplify=T)[,3],
#          strat_nAz=strat_Az*index) %>%
#   mutate(country=ifelse(bcr=="22"|bcr=="23"|bcr=="24", "Midwest (BCR 22, 23, 24)", "Other")) %>%
#   group_by(country, year, .draw) %>%
#   summarise(aggregate_N=sum(strat_nAz), aggregate_A=sum(strat_A)) %>%
#   ungroup()

# by country
samps4 <- Azs %>% 
  right_join(samps1) %>% 
  mutate(country=str_split(strat_name, "-", simplify=T)[,1],
         state=str_split(strat_name, "-", simplify=T)[,2],
         bcr=str_split(strat_name, "-", simplify=T)[,3],
         strat_nAz=strat_Az*index) %>% 
  group_by(country, year, .draw) %>% 
  summarise(aggregate_N=sum(strat_nAz), aggregate_A=sum(strat_A)) %>% 
  ungroup()

# plot the trajectories per country
pcntry <- samps4 %>% 
  ggplot(aes(x=year, y=aggregate_N/1000000)) +
  stat_lineribbon(.width=c(0.5, 0.8, 0.95), lwd=0.8) + 
  scale_fill_scico_d("Credible\ninterval", palette = "vik", 
                   begin=0.45, end=0.15) +
  facet_wrap(~country, scale="free") + theme_bw() +
  scale_x_continuous(breaks=seq(min(samps4$year), max(samps4$year), 4)) +
  labs(x="Year", y="Relative abundance")
plot_grid(pcntry, plot_grid(pall, p2008), ncol=1)

# poisson regression trends since 2008 per country
cbc_bcr_reg_trends_2008 <- samps4 %>% 
  filter(year>=2008) %>% 
  group_by(country, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(country) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))
cbc_bcr_reg_trends_2008

# poisson regression trends all years per country
cbc_bcr_reg_trends_1994 <- samps4 %>% 
  group_by(country, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(country) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))
cbc_bcr_reg_trends_1994
# ------------------------------------------------------------------------------


# aggregations and trends continent ----------------------------------------------
# aggregate the abundance indices across stratum and sample per continent
samps5 <- Azs %>% 
  right_join(samps1) %>% 
  mutate(country="North America",
         state=str_split(strat_name, "-", simplify=T)[,2],
         bcr=str_split(strat_name, "-", simplify=T)[,3],
         strat_nAz=strat_Az*index) %>% 
  group_by(country, year, .draw) %>% 
  summarise(aggregate_N=sum(strat_nAz), aggregate_A=sum(strat_A)) %>% 
  ungroup()

# plot the trajectories continent
pna <- samps5 %>% 
  ggplot(aes(x=year, y=aggregate_N/1000000)) +
  stat_lineribbon(.width=c(0.5, 0.8, 0.95), lwd=0.8) + 
  scale_fill_scico_d("Credible\ninterval", palette = "vik", 
                     begin=0.45, end=0.15) +
  facet_wrap(~country, scale="free") + theme_bw() +
  scale_x_continuous(breaks=seq(min(samps4$year), max(samps4$year), 4)) +
  labs(x="Year", y="Relative abundance")
plot_grid(pna)

# poisson regression trends since 2008 per country
cbc_cont_reg_trends_2008 <- samps5 %>% 
  filter(year>=2008, year<=2018) %>% 
  group_by(country, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(country) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))
cbc_cont_reg_trends_2008

# poisson regression trends all years per country
cbc_cont_reg_trends_1994 <- samps5 %>% 
  group_by(country, .draw) %>% 
  do(mod = glm(round(aggregate_N, 0) ~ year, data = ., family="poisson")) %>%
  mutate(slope = summary(mod)$coeff[2]) %>%
  select(-mod) %>% ungroup() %>% group_by(country) %>% 
  summarise(mean_trend=(exp(mean(slope))-1)*100, 
            lcl_trend=(exp(quantile(slope, probs=0.025))-1)*100, 
            ucl_trend=(exp(quantile(slope, probs=0.975))-1)*100) %>% 
  mutate(ci_range=ucl_trend-lcl_trend,
         diff_zero=ifelse(lcl_trend>0 | ucl_trend<0, 1, 0))
cbc_cont_reg_trends_1994
# ------------------------------------------------------------------------------


# # some other possibly models to try ------------------------------------------
# # bbs fixed intercept, fixed log-linear plus random year
# form <- count ~ 1 + firstyr + year_cent * factor(strat_idx1) +
#   f(obser_route_idx1, model="iid", constr=T) +
#   f(year_idx1, model="iid", group=strat_idx2, control.group=list(model="iid"))
# 
# # bbs random intercept, random log-linear plus random year
# form <- count ~ 1 + firstyr + year_cent +
#   f(obser_route_idx1, model="iid", constr=T) +
#   f(strat_idx1, model="iid", constr=T) +
#   f(strat_idx2, year_cent, model="iid", constr=T) +
#   f(year_idx1, model="iid", group=strat_idx3, control.group=list(model="iid"))
# 
# # bbs random intercept, first order year
# form <- count ~ 1 + firstyr + 
#   f(obser_route_idx1, model="iid", constr=T) + # observer-route effect
#   f(strat_idx1, model="iid", constr=T) + # intercept per stratum
#   f(year_idx1, model="rw1", constr=T, # rw (or 'ar1') year effect per stratum
#     group=strat_idx2, control.group=list(model="iid")) 
# 
# # cbc random intercept, random log-linear plus random year
# form <- count ~ 1 + log_hrs + year_cent +
#   f(circ_idx1, model="iid", constr=T) +
#   f(strat_idx1, log_hrs, model="iid", constr=T) +
#   f(strat_idx2, model="iid", constr=T) +
#   f(strat_idx3, year_cent, model="iid", constr=T) +
#   f(year_idx1, model="iid", constr=T, group=strat_idx4,
#     control.group=list(model="iid"))
# 
# # try em out
# test <- inla(form, data=d3,
#             family = "nbinomial", 
#             control.predictor=list(compute=T, link=1),
#             control.compute=list(config=T, waic=T, cpo=T,
#                                  control.gcpo=list(enable=TRUE)),
#             verbose=T, inla.mode="experimental",
#             control.inla = list(strategy = 'adaptive', int.strategy = 'eb'))
# 
# # view diagnostics
# dp9 <- ggplot(data.frame(x=test$cpo$pit[-pred_idx1])) +
#   geom_histogram(aes(x=x), col="white") + theme_bw()
# dp10 <- ggplot(data.frame(predicted=test$summary.fitted.values$mean, 
#        observed=d3$count), aes(x=observed, y=predicted)) + 
#   geom_point() + geom_abline(col="darkblue") + theme_bw()
# plot_grid(dp9, dp10)
# cor(test$summary.fitted.values$mean, d3$count, use="complete.obs")
# 
# # model scores
# get_bpic(test, pred_idx1)$`-2*sum(log(cpo_vec), na.rm=T)`
# bpic1$`-2*sum(log(cpo_vec), na.rm=T)`
# test$waic$waic
# waic1
# # ------------------------------------------------------------------------------

