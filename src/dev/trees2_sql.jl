--You can run the below set.hive commands if you like. They do not impact the running of the exposure file code
set hive.cli.print.header=true;
set hive.exec.dynamic.partition.mode=nonstrict; 
set hive.vectorized.execution.enabled = true; 
set hive.vectorized.execution.reduce.enabled = true; 
set hive.cbo.enable=true; 
set hive.compute.query.using.stats=true; 
set hive.stats.fetch.column.stats=true; 
set hive.stats.fetch.partition.stats=true; 
set hive.optimize.insert.dest.volume=true; 
set hive.auto.convert.join.noconditionaltask = true; 
set hive.auto.convert.join.noconditionaltask.size = 10000000;
set hive.plan.serialization.format=javaXML;
set hive.plan.serialization.format=kyro;

--Location where all of the tables for this project will be stored
create database if not exists jennie;

--Load in the placement and creative lookup information
--Data from these tables are joined to the placement_id and creative_id in the exposure table in order to create breaks used in the analysis
--Prior to loading in the data, confirm that there are no duplicate records in the tables and that each placement_id and creative_id only appears once
--The comScore lookup table for placement_id and creative_id is typically sent to the delivery team where they create the groups

drop table if exists jennie.placements6;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.placements6(
publisher string,
id string,
placement_nm string)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ',' LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/placements6';
LOAD DATA INPATH '/mapr/mapr04p/analytics0001/analytic_users/mddak/jo_placements6.csv' INTO TABLE jennie.placements6;

drop table if exists jennie.creatives6;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.creatives6(
id string,
creative_nm string)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ',' LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/creatives6';
LOAD DATA INPATH '/mapr/mapr04p/analytics0001/analytic_users/mddak/jo_creatives6.csv' INTO TABLE jennie.creatives6;

--We can now start the process of creating campaign specific exposure files
--The date range for the analysis comes from the Spec sheet.
--Some campaigns require just clientid, like SUN and CDW, because we collapse (group) the individual campaign_ids for them into just one for analysis
--Other campaigns like Hormel, have the same clientid, but require the exposure files to be created using the specific campaign_id

drop table jennie.jennie6_breaks_RF;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.jennie6_breaks_RF(
VENUE_DIM_KEY BIGINT,
TM_DIM_KEY_DAY BIGINT,
date1 STRING,
Trans_TIME STRING,
Trans_NUM DECIMAL(20,0),
STORE_NUM BIGINT,
CARD_NUM BIGINT,
UPC BIGINT,
ndc_upc_ind string,
QUANTITY double ,
WEIGHT double,
item_list_price double,
NET_PRICE double ,
CARD_DISCOUNT double ,
ad_discount double,
other_discount double,
ITEM_TYPE   string ,  
TENDER_TYPE   int ,
ad_event string,           
ad_version bigint,  
TM_DIM_KEY_WEEK BIGINT,
OPERATING_COMPANY STRING,
item_dim_key bigint,
UPC10 STRING,
UNITS INT,
CENTS INT,
BASELINE_UNITS DOUBLE ,
BASELINE_CENTS DOUBLE ,
FEATURE INT,
DISPLAY INT,
TOTL_PRC_REDUC INT,
registration_request_id BIGINT,
product_id_Jennieo6 INT,
EXPERIAN_ID BIGINT,
gross_imps bigint,
val_imps bigint,
placement_nm string,
creative_nm string,
publisher string)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ',' LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/jennie6_breaks_RF';

Insert Overwrite Table jennie.jennie6_breaks_RF
select distinct
'0' as VENUE_DIM_KEY,
'0' as TM_DIM_KEY_DAY,
a.date1,
'na' as Trans_TIME,
'0' as Trans_NUM,
'0' as STORE_NUM,
'0' as CARD_NUM,
'0' as UPC,
'na' as ndc_upc_ind,
'0' as QUANTITY,
'0' as WEIGHT,
'0' as item_list_price,
'0' as NET_PRICE,
'0' as CARD_DISCOUNT,
'0' as ad_discount,
'0' as other_discount,
'na' as ITEM_TYPE,
'0' as TENDER_TYPE,
'na' as ad_event,
'0' as ad_version,
a.iri_week as TM_DIM_KEY_WEEK,
b.retailer as OPERATING_COMPANY,
'0' as item_dim_key,
'na' as UPC10,
'0' as UNITS,
'0' as CENTS,
'0' as BASELINE_UNITS,
'0' as BASELINE_CENTS,
'0' as FEATURE,
'0' as DISPLAY,
'0' as TOTL_PRC_REDUC,
'0' as registration_request_id,
'0' as product_id_Jennieo6,
b.household_id as experian_id,
a.gross_imps ,
a.val_imps ,
c.placement_nm,
d.creative_nm, 
c.publisher
from daily.daily_unioned a
left join jennie.placements6 c on a.placement_id = c.id
left join jennie.creatives6 d on a.creative_id = d.id
inner JOIN WH_SUPPLMENTAL.HOUSEHOLD_DEPENDENT_ID_MAP b on a.lr_id = b.dependent_id 
where b.dependent_type = 'LVRAMP_ID' and b.household_type = 'EXPERIAN' and b.retailer = 'COMSCORE' and b.data_supplier = 'EXPERIAN' and b.current_state = 'MATCHED'
and a.clientid in('21884504') and iri_week < 1934
order by experian_id, date1;

set hive.cli.print.header=true; select * from jennie.jennie6_breaks limit 10;
select publisher, sum(gross_imps), sum(val_imps), count(distinct lr_id) from jennie.jennie6_breaks group by publisher;
select placement_nm, sum(gross_imps), sum(val_imps), count(distinct lr_id) from jennie.jennie6_breaks group by placement_nm;
select creative_nm, sum(gross_imps), sum(val_imps), count(distinct lr_id) from jennie.jennie6_breaks group by creative_nm;


-- create weekly impression table
drop table if exists jennie.jennie6_daily_roll_up;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.jennie6_daily_roll_up(
lr_id string,
date1 int,
placement_nm string,
creative_nm string,
publisher string,
gross_imps bigint)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ','  LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/jennie6_daily_roll_up';

Insert Overwrite Table jennie.jennie6_daily_roll_up
select
lr_id,
date1,
placement_nm,
creative_nm,
publisher,
sum(gross_imps) as gross_imps
from jennie.jennie6_breaks
group by lr_id, date1, placement_nm, creative_nm, publisher
order by lr_id, date1;

set hive.cli.print.header=true; select * from jennie.jennie6_daily_roll_up limit 10;

-- create individual columns for breaks
drop table if exists jennie.jennie6_daily_breaks;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.jennie6_daily_breaks(
lr_id string,
date1 int,
HTML5 int,
Static int,
Native_creative int,
Video int,
Content int,
Contextual int,
Direct int,
Prospecting int,
Retargeting int,
Behavorial int,
Native int,
Predictive int,
3rdparty int,
PMP int,
Buzzfeed int,
Hulu int,
PopSugarUS int,
RachaelRayMag int,
Shape int,
Blank int,
EatingWell int,
MeredithCorporation int,
RealSimple int,
TURNDSP int,
CookingLight int,
dataxu int,
Allrecipes int,
Amazon int,
MyFitnessPal int,
WomensHealth int,
Yummly int,
YouTube int,
gross_imps bigint)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ','  LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/jennie6_daily_breaks';

Insert Overwrite Table jennie.jennie6_daily_breaks
select
lr_id,
date1,
(CASE WHEN creative_nm = 'HTML5' THEN gross_imps ELSE 0 END) AS HTML5,
(CASE WHEN creative_nm = 'Static' THEN gross_imps ELSE 0 END) AS Static,
(CASE WHEN creative_nm = 'Native' THEN gross_imps ELSE 0 END) AS Native_creative,
(CASE WHEN creative_nm = 'Video' THEN gross_imps ELSE 0 END) AS Video,
(CASE WHEN placement_nm = 'Content' THEN gross_imps ELSE 0 END) AS Content,
(CASE WHEN placement_nm = 'Contextual' THEN gross_imps ELSE 0 END) AS Contextual,
(CASE WHEN placement_nm = 'Direct' THEN gross_imps ELSE 0 END) AS Direct,
(CASE WHEN placement_nm = 'Prospecting' THEN gross_imps ELSE 0 END) AS Prospecting,
(CASE WHEN placement_nm = 'Retargeting' THEN gross_imps ELSE 0 END) AS Retargeting,
(CASE WHEN placement_nm = 'Behavorial' THEN gross_imps ELSE 0 END) AS Behavorial,
(CASE WHEN placement_nm = 'Native' THEN gross_imps ELSE 0 END) AS Native,
(CASE WHEN placement_nm = 'Predictive' THEN gross_imps ELSE 0 END) AS Predictive,
(CASE WHEN placement_nm = '3rd party' THEN gross_imps ELSE 0 END) AS 3rdparty,
(CASE WHEN placement_nm = 'PMP' THEN gross_imps ELSE 0 END) AS PMP,
(CASE WHEN publisher = 'Buzzfeed' THEN gross_imps ELSE 0 END) AS Buzzfeed,
(CASE WHEN publisher = 'Hulu.com' THEN gross_imps ELSE 0 END) AS Hulu,
(CASE WHEN publisher = 'PopSugar US' THEN gross_imps ELSE 0 END) AS PopSugarUS,
(CASE WHEN publisher = 'Rachael Ray Mag' THEN gross_imps ELSE 0 END) AS RachaelRayMag,
(CASE WHEN publisher = 'Shape' THEN gross_imps ELSE 0 END) AS Shape,
(CASE WHEN publisher = 'Blank' THEN gross_imps ELSE 0 END) AS Blank,
(CASE WHEN publisher = 'Eating Well' THEN gross_imps ELSE 0 END) AS EatingWell,
(CASE WHEN publisher = 'Meredith Corporation' THEN gross_imps ELSE 0 END) AS MeredithCorporation,
(CASE WHEN publisher = 'Real Simple' THEN gross_imps ELSE 0 END) AS RealSimple,
(CASE WHEN publisher = 'TURN DSP' THEN gross_imps ELSE 0 END) AS TURNDSP,
(CASE WHEN publisher = 'Cooking Light' THEN gross_imps ELSE 0 END) AS CookingLight,
(CASE WHEN publisher = 'advertisers.dataxu.com' THEN gross_imps ELSE 0 END) AS dataxu,
(CASE WHEN publisher = 'Allrecipes.com' THEN gross_imps ELSE 0 END) AS Allrecipes,
(CASE WHEN publisher = 'Amazon.com' THEN gross_imps ELSE 0 END) AS Amazon,
(CASE WHEN publisher = 'MyFitnessPal' THEN gross_imps ELSE 0 END) AS MyFitnessPal,
(CASE WHEN publisher = 'Womens Health' THEN gross_imps ELSE 0 END) AS WomensHealth,
(CASE WHEN publisher = 'Yummly Inc' THEN gross_imps ELSE 0 END) AS Yummly,
(CASE WHEN publisher = 'YouTube' THEN gross_imps ELSE 0 END) AS YouTube,
gross_imps
from jennie.jennie6_daily_roll_up
order by lr_id, date1;

set hive.cli.print.header=true; select * from jennie.jennie6_weekly_breaks limit 5;


-- sum up to one row per week per ID
drop table if exists jennie.jennie6_daily_breaks_ru;
CREATE EXTERNAL TABLE IF NOT EXISTS jennie.jennie6_daily_breaks_ru(
VENUE_DIM_KEY BIGINT,
TM_DIM_KEY_DAY BIGINT,
date1 STRING,
Trans_TIME STRING,
Trans_NUM DECIMAL(20,0),
STORE_NUM BIGINT,
CARD_NUM BIGINT,
UPC BIGINT,
ndc_upc_ind string,
QUANTITY double ,
WEIGHT double,
item_list_price double,
NET_PRICE double ,
CARD_DISCOUNT double ,
ad_discount double,
other_discount double,
ITEM_TYPE   string ,  
TENDER_TYPE   int ,
ad_event string,           
ad_version bigint,  
TM_DIM_KEY_WEEK BIGINT,
OPERATING_COMPANY STRING,
item_dim_key bigint,
UPC10 STRING,
UNITS INT,
CENTS INT,
BASELINE_UNITS DOUBLE ,
BASELINE_CENTS DOUBLE ,
FEATURE INT,
DISPLAY INT,
TOTL_PRC_REDUC INT,
registration_request_id BIGINT,
product_id_Jennieo6 INT,
household_id BIGINT,
HTML5 int,
Static int,
Native_creative int,
Video int,
Content int,
Contextual int,
Direct int,
Prospecting int,
Retargeting int,
Behavorial int,
Native int,
Predictive int,
3rdparty int,
PMP int,
Buzzfeed int,
Hulu int,
PopSugarUS int,
RachaelRayMag int,
Shape int,
Blank int,
EatingWell int,
MeredithCorporation int,
RealSimple int,
TURNDSP int,
CookingLight int,
dataxu int,
Allrecipes int,
Amazon int,
MyFitnessPal int,
WomensHealth int,
Yummly int,
YouTube int,
gross_imps bigint)
ROW FORMAT DELIMITED FIELDS TERMINATED BY ','  LOCATION '/mapr/mapr04p/analytics0001/analytic_users/mddak/jennie/jennie6_daily_breaks_ru';

Insert Overwrite Table jennie.jennie6_daily_breaks_ru
select
'0' as VENUE_DIM_KEY,
'0' as TM_DIM_KEY_DAY,
a.date1,
'na' as Trans_TIME,
'0' as Trans_NUM,
'0' as STORE_NUM,
'0' as CARD_NUM,
'0' as UPC,
'na' as ndc_upc_ind,
'0' as QUANTITY,
'0' as WEIGHT,
'0' as item_list_price,
'0' as NET_PRICE,
'0' as CARD_DISCOUNT,
'0' as ad_discount,
'0' as other_discount,
'na' as ITEM_TYPE,
'0' as TENDER_TYPE,
'na' as ad_event,
'0' as ad_version,
'0' as TM_DIM_KEY_WEEK,
max(b.retailer) as OPERATING_COMPANY,
'0' as item_dim_key,
'na' as UPC10,
'0' as UNITS,
'0' as CENTS,
'0' as BASELINE_UNITS,
'0' as BASELINE_CENTS,
'0' as FEATURE,
'0' as DISPLAY,
'0' as TOTL_PRC_REDUC,
'0' as registration_request_id,
'0' as product_id_Jennieo6,
b.household_id,
sum(a.HTML5) as HTML5,
sum(a.Static) as Static,
sum(a.Native_creative) as Native_creative,
sum(a.Video) as Video,
sum(a.Content) as Content,
sum(a.Contextual) as Contextual,
sum(a.Direct) as Direct,
sum(a.Prospecting) as Prospecting,
sum(a.Retargeting) as Retargeting,
sum(a.Behavorial) as Behavorial,
sum(a.Native) as Native,
sum(a.Predictive) as Predictive,
sum(a.3rdparty) as 3rdparty,
sum(a.PMP) as PMP,
sum(a.Buzzfeed) as Buzzfeed,
sum(a.Hulu) as Hulu,
sum(a.PopSugarUS) as PopSugarUS,
sum(a.RachaelRayMag) as RachaelRayMag,
sum(a.Shape) as Shape,
sum(a.Blank) as Blank,
sum(a.EatingWell) as EatingWell,
sum(a.MeredithCorporation) as MeredithCorporation,
sum(a.RealSimple) as RealSimple,
sum(a.TURNDSP) as TURNDSP,
sum(a.CookingLight) as CookingLight,
sum(a.dataxu) as dataxu,
sum(a.Allrecipes) as Allrecipes,
sum(a.Amazon) as Amazon,
sum(a.MyFitnessPal) as MyFitnessPal,
sum(a.WomensHealth) as WomensHealth,
sum(a.Yummly) as Yummly,
sum(a.YouTube) as YouTube,
sum(a.gross_imps) as gross_imps
from jennie.jennie6_daily_breaks a
inner JOIN WH_SUPPLMENTAL.HOUSEHOLD_DEPENDENT_ID_MAP b on a.lr_id = b.dependent_id 
where b.dependent_type = 'LVRAMP_ID' and b.household_type = 'EXPERIAN' and b.retailer = 'COMSCORE' and b.data_supplier = 'EXPERIAN' and b.current_state = 'MATCHED'
group by household_id, date1
order by household_id, date1;

set hive.cli.print.header=true; select * from jennie.jennie6_weekly_breaks_ru limit 5;


-- QC CHECKS ARE BELOW
select count(*) as row_cnt, count(distinct household_id) as ids, sum(gross_imps) as gross from jennie.jennie6_weekly_breaks_ru;

-- impression checks
-- compare these against previous read if possible. If there is a previous read make sure there is an increase in IDs and impressions in the new read
-- for this specific read make sure that there is no increase in IDs or impressions between steps. Important step is to make sure that for the 3rd check that the row_cnt and ids are equal

select count(*) as row_cnt, count(distinct lr_id) as ids, sum(gross_imps) as gross, sum(val_imps) as viewable, min(date1) as st_date, max(date1) as end_date from daily.daily_unioned where clientid in('21884504') and iri_week < 1934;
select count(*) as row_cnt, count(distinct lr_id) as ids, sum(gross_imps) as gross, sum(val_imps) as viewable, min(iri_week) as st_wk, max(iri_week) as end_wk from jennie.jennie6_breaks;


                                                                    
# ----- GT TEST ----
                                                                    
                                                                    
select count(*) as row_cnt, 
       count(distinct lr_id) as ids, 
       sum(gross_imps) as gross, 
       sum(val_imps) as viewable, 
       min(date1) as st_date, 
       max(date1) as end_date 
from daily.daily_unioned 
where clientid in('21884504') 
and iri_week < 1934;                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
