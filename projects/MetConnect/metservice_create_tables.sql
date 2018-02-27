BEGIN TRANSACTION
CREATE TABLE MetConnect.dbo.RainFallPredictionSitesGrid (
       MetConnectID int NOT NULL,
       SiteString varchar(34),
       Office varchar(2),
       HydroTelPointNo int,
       TidedaID int,
       StartDate date,
       PRIMARY KEY (MetConnectID)
)

CREATE TABLE MetConnect.dbo.RainFallPredictionsGrid (
       MetConnectID int NOT NULL,
       PredictionDateTime datetime NOT NULL,
       ReadingDateTime datetime NOT NULL,
       HourlyRainfall decimal(18,2),
       PRIMARY KEY (MetConnectID,PredictionDateTime,ReadingDateTime),
       FOREIGN KEY (MetConnectID) REFERENCES MetConnect.dbo.RainFallPredictionSitesGrid(MetConnectID)
)

insert into RainFallPredictionSitesGrid values (1, '13_mile_bush', 'C', 1030, 313710, '2010-11-22'),(2, 'athurs_pass', 'C', 846, 219510, '2010-11-22'),(3, 'geraldine_forest', 'C', 3611, 410010, '2010-11-22'),(4, 'grasmere', 'C', 656, 311810, '2010-11-22'),(5, 'kakahu_bush', 'C', 3613, 411012, '2010-11-22'),(6, 'mt_francis', 'C', 3612, 419211, '2010-11-22'),(7, 'kanuka_hill', 'C', 654, 320010, '2010-11-22'),(8, 'ranger_stm', 'C', 652, 218810, '2010-11-22'),(9, 'woodbury', 'C', 3615, 410111, '2010-11-22'),(27, 'opihi_at_dobson', 'C', 3617, 309610, '2012-07-10'),(28, 'opihi_at_kimbell', 'C', 3618, 400710, '2012-07-10'),(29, 'opihi_at_rockwood', 'C', 3621, 400910, '2012-07-10'),(30, 'opuha_at_clayton', 'C', 3616, 309810, '2012-07-10'),(31, 'opuha_at_fox', 'C', 3614, 303811, '2012-07-10'),(32, 'rocky_gully', 'C', 3619, 401610, '2012-07-10'),(33, 'tengawai_at_mackenzie_pass', 'C', 3620, 403711, '2012-07-10'),(34, 'lees_valley', 'C', 2618, 321212, '2015-05-05'),(35, 'kanuka_hill', 'C', 654, 320010, '2015-05-05'),(36, 'okuku', 'C', 2621, 322410, '2015-05-05'),(37, 'Okuku_school', 'C', 2620, 322411, '2015-05-05'),(38, 'Townshend', 'C', 2623, 322110, '2015-05-05'),(39, 'Ashley_Gorge', 'C', 2619, 322211, '2015-05-05'),(40, 'Pig_flat', 'C', 2622, 321310, '2015-05-05'),(41, 'south_ashburton_at_boundary_creek', 'C', 5363, 315010, '2017-06-30'),(42, 'north_ashburton_at_cookies_hut', 'C', 5364, 314411, '2017-06-30'),(43, 'clyde_river_at_erewhon', 'C', 5365, 305810, '2017-06-30'),(44, 'pudding_hill_at_mt_hutt', 'C', 5367, 315510, '2017-06-30'),(45, 'south_ashburton_at_mt_somers', 'C', 5368, 316310, '2017-06-30'),(46, 'ashburton_at_turtons_saddle', 'C', 5369, 314412, '2017-06-30')


ROLLBACK
