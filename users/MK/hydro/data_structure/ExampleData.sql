BEGIN TRANSACTION

INSERT INTO Hydro.dbo.DataProviderMaster (DataProvider, Description) VALUES
	('ECan', 'Our place'), ('NIWA', 'Bigger place')
	
INSERT INTO Hydro.dbo.SiteMaster (Status, Hazards, Owner, DataProviderID, Description) VALUES
	('Opened', 'None', 'ECan', 1, 'A nice site'), ('Opened', 'Big rocks', 'Bob Smith', 2, 'Another nice site')

INSERT INTO Hydro.dbo.FeatureMaster (FeatureLongName, FeatureShortName, Description) VALUES
	('River', 'river', 'A flowing freshwater surface water body'), ('Aquifer', 'aq', 'A flowing freshwater groundwater water body'), 
	('Atmosphere', 'atmos', 'Water vapour above ground')

INSERT INTO Hydro.dbo.MtypeMaster (MtypeLongName, MtypeShortName, MtypeGroup, Units, Description) VALUES
	('Flow', 'flow', 'Quantity', 'meter**3/second', 'volume per time'), ('Water Level', 'wl', 'Quanitity', 'meter', 'The water level above an arbitrary datum'),
	('Nitrate Nitrogen', 'no3', 'Quality', 'mg/l', 'Nitrate nitrogen concentration')

INSERT INTO Hydro.dbo.MeasurementSourceMaster (MSourceLongName, MSourceShortName, Description) VALUES
	('Recorder', 'rec', 'Automatic recording device in the field'), ('Manual field', 'mfield', 'Manually measured in the field'), 
	('Manual Lab', 'mlab', 'Manually collected in the field and measured in the lab')	
	
INSERT INTO Hydro.dbo.QualityStateMaster (QualityStateLongName, QualityStateShortName, Description) VALUES
	('RAW', 'raw', 'Raw data directly from measurement source'), ('Quality controlled', 'qc', 'Quality controlled data')
	
INSERT INTO Hydro.dbo.FeatureAttrMaster (FeatureID, FeatureAttrName, Description) VALUES
	(1, 'Minimum Flow Site', 'True'), (2, 'Well Depth', 'Well depth in meters'),
	(2, 'Well Diameter', 'Well diameter in cm')
	
INSERT INTO Hydro.dbo.SiteFeatureAttr (SiteID, FeatureAttrID, Value) VALUES
	(1, 1, 'True'), (2, 2, '10'), (2, 3, '20')	

INSERT INTO Hydro.dbo.LoggingMethodMaster (LoggingMethodName, Description) VALUES
	('Discrete', 'A discrete instantaneous measurement'), ('Sum', 'An accumulated sum since the last measurement'), 
	('Mean', 'An aggregated mean since the last measurement')

INSERT INTO Hydro.dbo.TSParamMaster (TSParamName, Description) VALUES
	('Lab Name', 'The lab name'), ('Lab Method', 'The lab method')
	
INSERT INTO Hydro.dbo.SiteFeature (SiteID, FeatureID, SiteFeatureNumber, SiteFeatureName) VALUES
	(1, 1, '69607', 'Big Stream at gorge'), (1, 2, 'M35/8765', 'Mistery Well'), (2, 3, '123456', 'Sky limit')		

INSERT INTO Hydro.dbo.FeatureMtype (FeatureID, MtypeID, Units) VALUES
	(1, 2, 'meter'), (1, 1, 'meter**3/second'), (1, 3, 'mg/l'), (2, 2, 'meter'), (2, 3, 'mg/l')
	
INSERT INTO Hydro.dbo.FeatureMtypeSource (FeatureID, MtypeID, MeasurementSourceID, QualityStateID, DataSource) VALUES
	(1, 1, 1, 1, 'Hydrotel'), (1, 1, 1, 2, 'Hydstra'), (1, 1, 2, 2, 'Bgauging'), (1, 3, 2, 2, 'Hilltop'), (2, 2, 2, 2, 'Wells'),
	(2, 3, 2, 2, 'Hilltop')
	
INSERT INTO Hydro.dbo.TSDataSite (SiteFeatureID, MtypeID, MeasurementSourceID, QualityStateID, LoggingMethodID, TSGroup) VALUES
	(1, 1, 1, 1, 1, 'Numeric'), (1, 1, 1, 2, 3, 'Numeric'), (1, 1, 2, 2, 1, 'Numeric'), (1, 3, 2, 2, 1, 'Character'), (2, 2, 2, 2, 1, 'Numeric'),
	(2, 3, 2, 2, 1, 'Character')	
	
	
COMMIT