# Partners response

This repository contains the R scripts and accompanying data sets used to reproduce the analyses presented in the manuscript by Wojczulanis‑Jakubas et al., “Limited capacity for parental compensation
under experimental handicap in a small Arctic alcid, the Little Auk,” currently *submitted to Behavioral Ecology*. The repository includes three
scripts, each corresponding to a specific research question addressed in the manuscript. Each script is paired with the relevant dataset(s), as
indicated within the script documentation.

**Question 1**: Does the provisioning effort of a parent change when its
partner is experimentally burdened?\
*Script:* FA_partners_feeding.R *Data:* nfeeds_chr1_df.rds

**Question 2**: Do any behavioural effects persist after the burden is
removed? *Script:* FA_partners_carryover.R *Data:* nfeeds_chr2_df.rds,
nfeeds_chr12_df.rds

**Question 3:** Are there ecological limits that could constrain
compensatory behaviour? *Script:* FA_partners_simulation.R *Data:*
activity_controls_chr1.rds

---

Descriptions of the script contents are provided within each script.  
The variables included in the associated datasets are defined as presented
below.

**Data:** *nfeeds_chr1_df.rds* and *nfeeds_chr2_df.rds* These two sets
share the same structure and represent the same type of measurments -
number of feeding performed within 48 h time window. The only difference
between the two sets is that *nfeeds_chr1.rds* represent the parameter
recorded during early chick rearing period (when devices were on) while
the *nfeeds_chr1.rds* represents the parameter during mid chick rearing
(when devices were off).

*season* - year of the study

*nest* - nest/pair identity

*sx* - sex of an individual

*nest_device_stat* - status of the nest in terms of experimental
treatment (control, devise: acc - accelerometer, gps - gps)

*partner_device_stat* - status of the pair member in respect to the
experimental treatment (control, deploy, partner)

*nfeeds* - number of feedings (i.e. foraging trips terminated with
feeding) performed by an individuals during 48h

---

**Data:** *nfeeds_chr12_df.rds* Same data as *nfeeds_chr1_df.rds* and
*nfeeds_chr2_df.rds*, just binded by rows to compare the birds
foraging/feeding performance between the two chick rearing periods
(early and mid; devices on and off)

*session* - recording session/chick rearing period: chick rearing1 -
early chick rearing (when devices were on) chick rearing2 - mid chick
rearing (when devices were off)

*season* - year of the study

*nest* - nest/pair identity

*sx* - sex of an individual

*nest_device_stat* - status of the nest in terms of experimental
treatment (control, devise: acc - accelerometer, gps - gps)

*partner_device_stat* - status of the pair member in respect to the
experimental treatment (control, deploy, partner)

*nfeeds* - number of feedings (i.e. foraging trips terminated with
feeding) performed by an individuals during 48h

---

**Data:** *activity_controls_chr1.rds* Data of complete time intervals
for the distinguished activities (see definition of those in the
manuscript), registered during early chick rearing period in control
birds.

*season* - year of the study

*ringno* - bird identity\

*behav_blocks* - distinguished activity (4 levels: foraging_flight,
latency, nest, colony - see definitions in the manuscript)
*tot_dur_min* - duration of the full time intervals for the
distinguished activities (in minutes)
